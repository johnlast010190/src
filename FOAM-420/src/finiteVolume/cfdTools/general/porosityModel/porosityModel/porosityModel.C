/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvc.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "coordinate/systems/cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityModel, 0);
    defineRunTimeSelectionTable(porosityModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::porosityModel::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorInFunction
            << "Negative resistances are invalid, resistance = " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


Foam::label Foam::porosityModel::fieldIndex(const label i) const
{
    label index = 0;
    if (!csys().uniform())
    {
        index = i;
    }
    return index;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel::porosityModel
(
    const word& name,
    const word& modelType,
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    obr_(obr),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(true),
    porousRelVelocity_(dict.lookupOrDefault<bool>("relativeVelocity", false)),
    zoneName_(cellZoneName),
    cellZoneIDs_(),
    coorFramePtr_(nullptr)
{
    if (zoneName_ == word::null)
    {
        dict.readIfPresent("active", active_);
        dict_.lookup("cellZone") >> zoneName_;
    }

    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    Info<< "    creating porous zone: " << zoneName_ << endl;

    bool foundZone = !cellZoneIDs_.empty();
    reduceToMaster(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorInFunction
            << "cannot find porous cellZone " << zoneName_
            << exit(FatalError);
    }

    coordinateSystem::errorCoordinateSystem(coeffs_, porosityModel::typeName);
    coorFramePtr_ =
        coeffs_.found("referenceFrame")
      ? coordinateFrame::lookupNew(mesh, coeffs_)
      : coordinateFrame::globalFrame(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModel::~porosityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::porosityModel::transformModelData()
{
    if (!mesh_.upToDatePoints(*this))
    {
        calcTransformModelData();

        // set model up-to-date wrt points
        mesh_.setUpToDatePoints(*this);
    }
}


Foam::tmp<Foam::vectorField> Foam::porosityModel::porosityModel::force
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    transformModelData();

    tmp<vectorField> tforce(new vectorField(U.size(), Zero));

    if (!cellZoneIDs_.empty())
    {
        this->calcForce(U, rho, mu, tforce.ref());
    }

    return tforce;
}


void Foam::porosityModel::addResistance(fvVectorMatrix& UEqn)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn);
}


void Foam::porosityModel::addResistance(fvBlockMatrix<vector>& UEqn)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn);
}


void Foam::porosityModel::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}


void Foam::porosityModel::addResistance
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}



void Foam::porosityModel::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, AU);

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}

void Foam::porosityModel::addAdjointResistance
(
    fvVectorMatrix& UaEqn,
    const volVectorField& Uprimal
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->adjointCorrect(UaEqn, Uprimal);
}


void Foam::porosityModel::addAdjointResistance
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->adjointCorrect(UaEqn, rho, mu, Uprimal);
}


void Foam::porosityModel::addAdjointResistance
(
    fvBlockMatrix<vector>& UaEqn,
    const volVectorField& U
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->adjointCorrect(UaEqn, U);
}



void Foam::porosityModel::addAdjointResistance
(
    fvBlockMatrix<vector>& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& U
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->adjointCorrect(UaEqn, rho, mu, U);
}


Foam::labelList Foam::porosityModel::zoneBoundaryFaces() const
{

    if (cellZoneIDs_.empty())
    {
        return labelList(0);
    }

    DynamicList<label> zoneBoundFaces(100);


    boolList zoneCells(mesh_.nCells(), false);
    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            zoneCells[cells[i]] = true;
        }
    }

    const labelList& nei = mesh_.neighbour();
    const labelList& own = mesh_.faceOwner();

    forAll(nei, i)
    {
        if (zoneCells[own[i]] && !zoneCells[nei[i]])
    {
            zoneBoundFaces.append(i);
    }
    else if (zoneCells[nei[i]] && !zoneCells[own[i]])
    {
            zoneBoundFaces.append(i);
    }
    }

    // check on processor patches
    boolList neiZoneID(mesh_.nFaces()-mesh_.nInternalFaces(),false);

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        if (zoneCells[own[faceI]])
        {
            neiZoneID[faceI-mesh_.nInternalFaces()] =  true;
        }
    }
    boolList swappedZoneID(neiZoneID);

    syncTools::swapBoundaryFaceList(mesh_, swappedZoneID);

    const polyPatchList& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            forAll(procPatch, i)
            {
                label meshFaceI = procPatch.start() + i;
                label bFaceI =  meshFaceI - mesh_.nInternalFaces();

                if
                (
                   !swappedZoneID[bFaceI]
                 && neiZoneID[bFaceI]
                )
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
        else
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            forAll(patch, i)
            {
                label meshFaceI = patch.start() + i;
                label bFaceI = meshFaceI - mesh_.nInternalFaces();
                if (neiZoneID[bFaceI])
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
    }

    return labelList(zoneBoundFaces, true);
}


bool Foam::porosityModel::writeData(Ostream& os) const
{
    return true;
}


bool Foam::porosityModel::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    coeffs_ = dict.optionalSubDict(type() + "Coeffs");

    dict.lookup("cellZone") >> zoneName_;
    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    return true;
}


// ************************************************************************* //
