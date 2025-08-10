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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "thrusterActuationDiskSource.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::fv::
        thrusterActuationDiskSource::fanModelType,
        3
    >::names[] =
    {
        "averageVelocity",
        "localVelocity",
        "flowRate"
    };

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(thrusterActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        thrusterActuationDiskSource,
        dictionary
    );

//
const NamedEnum
<
    Foam::fv::thrusterActuationDiskSource::fanModelType,
    3
> Foam::fv::thrusterActuationDiskSource::fanModelTypeNames_;

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::thrusterActuationDiskSource::thrusterActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    fanModel_(fanModelTypeNames_.read(coeffs_.lookup("fanMode"))),
    load_(cells_.size(), Zero),
    jumpV_(cells_.size(), vector::zero),
    f_(coeffs_.lookup("f")),
    lowerBound_(coeffs_.lookupOrDefault("lowerBound", 0.0)),
    upperBound_(coeffs_.lookupOrDefault("upperBound", 1E+06)),
    coorFramePtr_(coordinateFrame::lookupNew(mesh_, coeffs_)),
    loadProfile_(Function1<scalar>::New("loadProfile", coeffs_)),
    alphaProfile_(Function1<scalar>::New("alphaProfile", coeffs_)),
    zoneBoundaryFaces_()
{
    Info<< "    - creating thruster actuation disk zone: " << name_ << endl;

    fieldNames_ = coeffs_.lookup<wordList>("fieldNames");
    applied_.setSize(fieldNames_.size(), false);

    // Incompressible density
    rhoRefFound_ = coeffs_.found("rhoRef");
    if (rhoRefFound_)
    {
        coeffs_.lookup("rhoRef") >> rhoRef_;
    }
    else
    {
        rhoRef_ = 1.0;
    }

    if (mag(coorFramePtr_->axis()) < SMALL)
    {
        FatalErrorInFunction
            << "Badly defined axis: zero magnitude: " << coorFramePtr_->axis()
            << abort(FatalError);
    }
    coorFramePtr_->axis() /= mag(coorFramePtr_->axis());

    scalarField rOverR(cells_.size());
    vectorField tVel(cells_.size());

    PackedBoolList sourcePts(mesh_.nPoints(), 0);

    forAll(cells_, i)
    {
        label celli = cells_[i];
        const point cc = mesh_.cellCentres()[celli];
        vector v(cc - coorFramePtr_->CofR());
        scalar parallel = (v & coorFramePtr_->axis());

        //remove parallel component
        v -= (parallel*coorFramePtr_->axis());

        rOverR[i] = mag(v);
        tVel[i] = v ^ coorFramePtr_->axis();

        const labelList& cellPoints = mesh_.cellPoints()[celli];
        forAll(cellPoints, cPtI)
        {
            sourcePts.set(cellPoints[cPtI], 1);
        }
    }

    //Calculate source thickness
    scalar minThickness = GREAT;
    scalar maxThickness = -GREAT;
    forAll(mesh_.points(), pointi)
    {
        if (sourcePts.get(pointi) == 1)
        {
            vector v(mesh_.points()[pointi] - coorFramePtr_->CofR());
            scalar parallel = (v & coorFramePtr_->axis());
            minThickness = min(minThickness, parallel);
            maxThickness = max(maxThickness, parallel);
        }
    }
    reduce(minThickness, minOp<scalar>());
    reduce(maxThickness, maxOp<scalar>());

    scalar thickness = mag(maxThickness - minThickness);
    Info<< "Thruster thickness : " << thickness << endl;

    diskArea_ = V_/thickness;
    Info<< "Disk Area : " << diskArea_ << endl;

    const scalar maxR = gMax(rOverR);
    rOverR /= maxR;

    load_ = loadProfile_->value(rOverR)/thickness;

    const scalar fs = mag(coorFramePtr_->Omega());

    jumpV_ =
        fs*tVel*Foam::sin
        (
            alphaProfile_->value(rOverR)
           /180.*Foam::constant::mathematical::pi
        )/thickness;

    zoneBoundaryFaces_ = zoneBoundaryFaces();

    compressible_ = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fv::thrusterActuationDiskSource::zoneBoundaryFaces() const
{
    DynamicList<label> zoneBoundFaces(0);

    boolList zoneCells(mesh_.nCells(), false);

    forAll(cells_, i)
    {
        zoneCells[cells_[i]] = true;
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
    boolList neiZoneID(mesh_.nFaces() - mesh_.nInternalFaces(), false);

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        if (zoneCells[own[facei]])
        {
            neiZoneID[facei - mesh_.nInternalFaces()] =  true;
        }
    }

    boolList swappedZoneID(neiZoneID);

    syncTools::swapBoundaryFaceList(mesh_, swappedZoneID);

    const polyPatchList& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            forAll(procPatch, i)
            {
                const label facei = procPatch.start() + i;
                const label bFaceI = facei - mesh_.nInternalFaces();

                if
                (
                    !swappedZoneID[bFaceI] && neiZoneID[bFaceI]
                )
                {
                    zoneBoundFaces.append(facei);
                }
            }
        }
        else
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchi];
            forAll(patch, i)
            {
                label facei = patch.start() + i;
                label bFaceI =  facei - mesh_.nInternalFaces();
                if (neiZoneID[bFaceI])
                {
                    zoneBoundFaces.append(facei);
                }
            }
        }
    }

    return labelList(zoneBoundFaces, true);
}


Foam::scalar Foam::fv::thrusterActuationDiskSource::zoneFlux
(
    const fvMesh& mesh,
    const labelList& faces
) const
{
    scalar fsum = 0;
    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    forAll(faces, i)
    {
        const label facei = faces[i];
        if (facei < mesh.nInternalFaces())
        {
            fsum += mag(phi[facei]);
        }
        else
        {
            const label patchi = mesh.boundaryMesh().whichPatch(facei);
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            if (!isA<emptyPolyPatch>(patch))
            {
                label bfI = facei - patch.start();
                fsum += mag(phi.boundaryField()[patchi][bfI]);
            }
        }
    }
    reduce(fsum, sumOp<scalar>());
    return fsum;
}


void Foam::fv::thrusterActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        compressible_ = true;
    }

    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        if (compressible_)
        {
            addThrusterActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                obr_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            if (!rhoRefFound_)
            {
                FatalErrorInFunction
                    << "rhoRef must be specified for incompressible flow."
                    << exit(FatalError);
            }
            addThrusterActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                geometricOneField(),
                U
            );
        }
    }
}


void Foam::fv::thrusterActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        compressible_ = true;
    }

    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        if (compressible_)
        {
            addThrusterActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                rho,
                U
            );
        }
        else
        {
            if (!rhoRefFound_)
            {
                FatalErrorInFunction
                    << "rhoRef must be specified for incompressible flow."
                    << exit(FatalError);
            }
            addThrusterActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                geometricOneField(),
                U
            );
        }
    }
}


void Foam::fv::thrusterActuationDiskSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::thrusterActuationDiskSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
