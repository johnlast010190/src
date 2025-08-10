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

#include "finiteVolume/gradSchemes/limitedGradSchemes/cellMDLimitedGrad/cellMDLimitedGrad.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvGradScheme(cellMDLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::cellMDLimitedGrad<Foam::scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volVectorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalarField maxVsf(vsf.primitiveField());
    scalarField minVsf(vsf.primitiveField());

    surfaceScalarField mvsfOwn(fvc::ownerField(vsf));
    surfaceScalarField mvsfNei(fvc::neighbourField(vsf));
    fvc::applyFaceMaskTo(mvsfOwn);
    fvc::applyFaceMaskTo(mvsfNei);
    // Add and sync the masked indirect and internal faces in both fields
    fvc::sumIndirectAndInternalFaces(mvsfOwn, mvsfNei);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        maxVsf[own] = max(maxVsf[own], mvsfNei[facei]);
        minVsf[own] = min(minVsf[own], mvsfNei[facei]);

        maxVsf[nei] = max(maxVsf[nei], mvsfOwn[facei]);
        minVsf[nei] = min(minVsf[nei], mvsfOwn[facei]);
    }


    const volScalarField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        const scalarField& pmvsfNei(mvsfNei.boundaryField()[patchi]);
        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            maxVsf[own] = max(maxVsf[own], pmvsfNei[pFacei]);
            minVsf[own] = min(minVsf[own], pmvsfNei[pFacei]);
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const scalarField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            g[own],
            maxVsf[own],
            minVsf[own],
            Cf[facei] - C[own]
        );

        // neighbour side
        limitFace
        (
            g[nei],
            maxVsf[nei],
            minVsf[nei],
            Cf[facei] - C[nei]
        );
    }


    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                g[own],
                maxVsf[own],
                minVsf[own],
                pCf[pFacei] - C[own]
            );
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::cellMDLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    vectorField maxVsf(vsf.primitiveField());
    vectorField minVsf(vsf.primitiveField());

    surfaceVectorField mvsfOwn(fvc::ownerField(vsf));
    surfaceVectorField mvsfNei(fvc::neighbourField(vsf));
    fvc::applyFaceMaskTo(mvsfOwn);
    fvc::applyFaceMaskTo(mvsfNei);
    // Add and sync the masked indirect and internal faces in both fields
    fvc::sumIndirectAndInternalFaces(mvsfOwn, mvsfNei);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        maxVsf[own] = max(maxVsf[own], mvsfNei[facei]);
        minVsf[own] = min(minVsf[own], mvsfNei[facei]);

        maxVsf[nei] = max(maxVsf[nei], mvsfOwn[facei]);
        minVsf[nei] = min(minVsf[nei], mvsfOwn[facei]);
    }


    const volVectorField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        const vectorField& pmvsfNei(mvsfNei.boundaryField()[patchi]);
        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            maxVsf[own] = max(maxVsf[own], pmvsfNei[pFacei]);
            minVsf[own] = min(minVsf[own], pmvsfNei[pFacei]);
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const vectorField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            g[own],
            maxVsf[own],
            minVsf[own],
            Cf[facei] - C[own]
        );

        // neighbour side
        limitFace
        (
            g[nei],
            maxVsf[nei],
            minVsf[nei],
            Cf[facei] - C[nei]
        );
    }


    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                g[own],
                maxVsf[own],
                minVsf[own],
                pCf[pFacei] - C[own]
            );
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


// ************************************************************************* //
