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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "singleLayerRegion/singleLayerRegion.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(singleLayerRegion, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::singleLayerRegion::constructMeshObjects()
{
    // construct patch normal vectors
    nHatPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "nHat",
                time_.timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                NO_WRITE
            ),
            regionMesh(),
            dimensionedVector("zero", dimless, Zero),
            zeroGradientFvPatchField<vector>::typeName
        )
    );

    // construct patch areas
    magSfPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "magSf",
                time_.timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                NO_WRITE
            ),
            regionMesh(),
            dimensionedScalar("zero", dimArea, 0.0),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
}


void Foam::regionModels::singleLayerRegion::initialise()
{
    if (debug)
    {
        Pout<< "singleLayerRegion::initialise()" << endl;
    }

    label nBoundaryFaces = 0;
    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();
    volVectorField& nHat = nHatPtr_();
    volScalarField& magSf = magSfPtr_();
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = rbm[patchi];
        const labelList& fCells = pp.faceCells();

        nBoundaryFaces += fCells.size();

        UIndirectList<vector>(nHat, fCells) = pp.faceNormals();
        UIndirectList<scalar>(magSf, fCells) = pp.magFaceAreas();
    }
    nHat.correctBoundaryConditions();
    magSf.correctBoundaryConditions();

    if (nBoundaryFaces != regionMesh().nCells())
    {
        FatalErrorInFunction
            << "Number of primary region coupled boundary faces not equal to "
            << "the number of cells in the local region" << nl << nl
            << "Number of cells = " << regionMesh().nCells() << nl
            << "Boundary faces  = " << nBoundaryFaces << nl
            << abort(FatalError);
    }

    scalarField passiveMagSf(magSf.size(), 0.0);
    passivePatchIDs_.setSize(intCoupledPatchIDs_.size(), -1);
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppIntCoupled = rbm[patchi];
        if (ppIntCoupled.size() > 0)
        {
            label cellId = rbm[patchi].faceCells()[0];
            const cell& cFaces = regionMesh().cells()[cellId];

            label facei = ppIntCoupled.start();
            label faceO = cFaces.opposingFaceLabel(facei, regionMesh().faces());

            label passivePatchi = rbm.whichPatch(faceO);
            passivePatchIDs_[i] = passivePatchi;
            const polyPatch& ppPassive = rbm[passivePatchi];
            UIndirectList<scalar>(passiveMagSf, ppPassive.faceCells()) =
                ppPassive.magFaceAreas();
        }
    }

    Pstream::listCombineReduce(passivePatchIDs_, maxOp<label>());

    magSf.field() = 0.5*(magSf + passiveMagSf);
    magSf.correctBoundaryConditions();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::singleLayerRegion::read()
{
    return regionModel::read();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::singleLayerRegion
(
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel(mesh, regionType),
    nHatPtr_(nullptr),
    magSfPtr_(nullptr),
    passivePatchIDs_()
{}


Foam::regionModels::singleLayerRegion::singleLayerRegion
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    nHatPtr_(nullptr),
    magSfPtr_(nullptr),
    passivePatchIDs_()
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::~singleLayerRegion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::regionModels::singleLayerRegion::nHat() const
{
    if (!nHatPtr_.valid())
    {
        FatalErrorInFunction
            << "Region patch normal vectors not available"
            << abort(FatalError);
    }

    return nHatPtr_();
}


const Foam::volScalarField& Foam::regionModels::singleLayerRegion::magSf() const
{
    if (!magSfPtr_.valid())
    {
        FatalErrorInFunction
            << "Region patch areas not available"
            << abort(FatalError);
    }

    return magSfPtr_();
}


const Foam::labelList&
Foam::regionModels::singleLayerRegion::passivePatchIDs() const
{
    return passivePatchIDs_;
}


// ************************************************************************* //
