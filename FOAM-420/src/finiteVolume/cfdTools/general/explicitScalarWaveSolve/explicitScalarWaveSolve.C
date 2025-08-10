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
    (c) Esi Ltd

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/explicitScalarWaveSolve/explicitScalarWaveSolve.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvc.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveSolverCoeffs::waveSolverCoeffs
(
    const fvMesh& mesh,
    const vector& direction
)
:
    labelList(mesh.nCells(), -1),
    A_(mesh.nCells(), 0),
    HCoeffs_
    (
        mesh.nCells(),
        List<Tuple2<scalar, label>>(0, Tuple2<scalar, label>(0.0, -1))
    ),
    Jppos_
    (
        mesh.boundary(),
        mesh.V(),           // Dummy internal field,
        calculatedFvPatchScalarField::typeName
    )
{
    // base flux field
    surfaceScalarField Ji(direction & mesh.Sf());

    // face flux flipping
    const labelList& allOwner = mesh.faceOwner();
    const labelList& allNeigh = mesh.faceNeighbour();

    // calculate optimal solution order
    {
        SortableList<scalar> dist
        (
            scalarField(direction & mesh.C().internalField())
        );

        this->labelList::operator=(labelList(dist.indices()));


        labelList inverseOrder(dist.size(), -1);

        forAll(inverseOrder, ci)
        {
            inverseOrder[this->operator[](ci)] = ci;
        }

        label swapped = 0;
        label swapped_0 = -1;
        label nSwapsEqual = 0;
        do
        {
            //check for infinite loops
            if (swapped_0 == swapped)
            {
                nSwapsEqual++;
            }
            swapped_0 = swapped;
            swapped = 0;

            //try bubble-like sort based on face fluxes
            forAll(allNeigh, fi)
            {
                label neiPos = inverseOrder[allNeigh[fi]];
                label ownPos = inverseOrder[allOwner[fi]];

                if (Ji[fi] != 0)
                {
                    if (sign(ownPos - neiPos) == sign(Ji[fi]))
                    {
                        //swap!
                        swapped++;

                        label neiPosCell = this->operator[](neiPos);
                        this->operator[](neiPos) = this->operator[](ownPos);
                        this->operator[](ownPos) = neiPosCell;
                        inverseOrder[allNeigh[fi]] = ownPos;
                        inverseOrder[allOwner[fi]] = neiPos;
                    }
                }
            }

        } while (swapped > 0 && nSwapsEqual < 10);

    }

    // now build coefficients and boundary field
    forAll(*this, ici)
    {
        label cellI = this->operator[](ici);

        const labelList& cellFaces(mesh.cells()[cellI]);

        DynamicList<scalar> dwFluxes(cellFaces.size());
        DynamicList<label> dwCells(cellFaces.size());

        forAll(cellFaces, cfI)
        {
            scalar Jif = 0;
            const label gFaceI = cellFaces[cfI];

            if (gFaceI < mesh.nInternalFaces())
            {
                Jif = Ji[gFaceI];
            }
            else
            {
                const label patchI
                    = mesh.boundaryMesh().whichPatch(gFaceI);
                const label startFace
                    = mesh.boundaryMesh()[patchI].start();
                Jif = Ji.boundaryField()[patchI][gFaceI-startFace];
            }


            bool isOwner (allOwner[gFaceI] == cellI);

            if ((Jif < 0) == isOwner)
            {
                if (gFaceI < mesh.nInternalFaces())
                {
                    dwFluxes.append(mag(Ji[gFaceI]));

                    if (isOwner)
                    {
                        dwCells.append(allNeigh[gFaceI]);
                    }
                    else
                    {
                        dwCells.append(allOwner[gFaceI]);
                    }
                }
                else
                {
                    const label patchI
                        = mesh.boundaryMesh().whichPatch(gFaceI);
                    const label startFace
                        = mesh.boundaryMesh()[patchI].start();
                    Jppos_[patchI][gFaceI-startFace] = mag(Jif);
                }
            }
            else
            {
                A_[ici] += mag(Jif);
            }

        }

        // assign to compact format
        dwFluxes.shrink();
        dwCells.shrink();

        HCoeffs_[ici].setSize(dwFluxes.size());

        forAll(HCoeffs_[ici], coi)
        {
            HCoeffs_[ici][coi].first() = dwFluxes[coi];
            HCoeffs_[ici][coi].second() = dwCells[coi];
        }

    }

}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


scalar explicitScalarWaveSolve
(
    volScalarField& vsf,
    surfaceScalarField& ssf,
    const surfaceScalarField& flux,
    const labelList& order
)
{
    const fvMesh& mesh = vsf.mesh();

    // find maximum intensity on boundary
    // internal values cannot be higher than this
    scalar vsfmax = 0;

    forAll(vsf.boundaryField(), bI)
    {
        vsfmax = max(vsfmax, max(vsf.boundaryField()[bI]));
    }
    scalar maxEqRes = 0;

    const labelList& allOwner = mesh.faceOwner();

    forAll(order, soI)
    {
        const label cI = order[soI];

        scalar& Ic = vsf[cI];

        const labelList& faces = mesh.cells()[cI];

        scalar lhs = 0;
        scalar rhs = 0;

        forAll(faces, fI)
        {
            const label gFaceI = faces[fI];

            //get flux and value for this face
            scalar Jif = 0;
            scalar If = 0;

            if (gFaceI < mesh.nInternalFaces())
            {
                Jif = flux[gFaceI];
                If = ssf[gFaceI];
            }
            else
            {
                const label patchI
                    = mesh.boundaryMesh().whichPatch(gFaceI);
                const label startFace
                    = mesh.boundaryMesh()[patchI].start();

                if (!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI]))
                {
                    Jif = flux.boundaryField()[patchI][gFaceI-startFace];
                    If = ssf.boundaryField()[patchI][gFaceI-startFace];
                }
            }

            if (allOwner[gFaceI] == cI)
            {
                Jif *= -1;
            }

            if (Jif > 0)
            {
                rhs += Jif * If;
            }
            else
            {
                lhs -= Jif;
            }
        }

        scalar resErr = mag(rhs - lhs*Ic)/max(rhs, VSMALL);
        Ic = rhs / lhs;
        Ic = min(Ic, vsfmax);
        Ic = max(Ic, 0.0);

        //update internal face values using UD and latest cell value
        forAll(faces, fI)
        {
            const label gFaceI = faces[fI];
            //get flux and value for this face
            if (gFaceI < mesh.nInternalFaces())
            {
                const scalar Jif = flux[gFaceI];
                if (allOwner[gFaceI] == cI && Jif > 0)
                {
                    ssf[gFaceI] = Ic;
                }
                else if (Jif <= 0)
                {
                    ssf[gFaceI] = Ic;
                }
            }
        }

        maxEqRes = max(resErr, maxEqRes);
    }

    reduce(maxEqRes, maxOp<scalar>());

    Info<< "Explicit Scalar Wave: Solving for "
         << vsf.name()
         << ", Maximum initial residual = " << maxEqRes << endl;

    return maxEqRes;
}

scalar explicitScalarWaveSolve
(
    volScalarField& vsf,
    const waveSolverCoeffs& wsc
)
{
    //intermediate solution field for boundary flux contributions
    scalarField tvsf(vsf.internalField().size(), 0.0);

    // find maximum intensity on boundary
    // internal values cannot be higher than this
    scalar vsfmax = 0;

    forAll(vsf.boundaryField(), bI)
    {
        vsfmax = max(vsfmax, max(vsf.boundaryField()[bI]));
    }
    scalar maxEqRes = 0;


    //add boundary contributions first
    forAll(vsf.boundaryField(), pi)
    {
        const labelUList& patchLabels
            (vsf.boundaryField()[pi].patch().faceCells());

        forAll(patchLabels, pfi)
        {

            tvsf[patchLabels[pfi]]
                += vsf.boundaryField()[pi][pfi]
                * wsc.Jppos_[pi][pfi];
        }
    }

    //main cell value calculation loop
    forAll(wsc, ici)
    {

        label ci = wsc[ici];


        forAll(wsc.HCoeffs_[ici], dfi)
        {
            tvsf[ci]
                += wsc.HCoeffs_[ici][dfi].first()
                * vsf[wsc.HCoeffs_[ici][dfi].second()];
        }

        scalar resErr
            = mag(tvsf[ci] - wsc.A_[ici]*vsf[ci])
            /max(tvsf[ci], VSMALL);


        vsf[ci] = tvsf[ci] / wsc.A_[ici];

        // bound
        vsf[ci] = max(0, vsf[ci]);
        vsf[ci] = min(vsf[ci], vsfmax);

        //max err
        maxEqRes = max(resErr, maxEqRes);

    }

    Info<< "Explicit Scalar Wave: Solving for "
         << vsf.name()
         << ", Maximum initial residual = " << maxEqRes << endl;

    return maxEqRes;
}


// ************************************************************************* //
} //end namespace FOAM

// ************************************************************************* //
