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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/limitedGradSchemes/faceLimitedGrad/faceLimitedGrad.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvGradScheme(faceLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::faceLimitedGrad<Foam::scalar>::calcGrad
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

    // create limiter
    scalarField limiter(vsf.primitiveField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

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

        scalar maxFace = max(mvsfOwn[facei], mvsfNei[facei]);
        scalar minFace = min(mvsfOwn[facei], mvsfNei[facei]);
        scalar maxMinFace = rk*(maxFace - minFace);
        maxFace += maxMinFace;
        minFace -= maxMinFace;

        // owner side
        limitFace
        (
            limiter[own],
            maxFace - vsf[own], minFace - vsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxFace - vsf[nei], minFace - vsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    const volScalarField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        const scalarField& pmvsfOwn(mvsfOwn.boundaryField()[patchi]);
        const scalarField& pmvsfNei(mvsfNei.boundaryField()[patchi]);

        if (psf.coupled() || psf.fixesValue())
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                scalar maxFace = max(pmvsfOwn[pFacei], pmvsfNei[pFacei]);
                scalar minFace = min(pmvsfOwn[pFacei], pmvsfNei[pFacei]);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsf[own], minFace - vsf[own],
                    (pCf[pFacei] - C[own]) & g[own]
                );
            }
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.primitiveFieldRef() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::faceLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vvf,
    const word& name
) const
{
    const fvMesh& mesh = vvf.mesh();

    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vvf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    // create limiter
    scalarField limiter(vvf.primitiveField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    surfaceVectorField vvfOwn(fvc::ownerField(vvf));
    surfaceVectorField vvfNei(fvc::neighbourField(vvf));
    tmp<surfaceVectorField> tmvvfOwn = fvc::applyFaceMask(vvfOwn);
    tmp<surfaceVectorField> tmvvfNei = fvc::applyFaceMask(vvfNei);
    surfaceVectorField& mvvfOwn = const_cast<surfaceVectorField&>(tmvvfOwn());
    surfaceVectorField& mvvfNei = const_cast<surfaceVectorField&>(tmvvfNei());
    // Add and sync the masked indirect and internal faces in both fields
    fvc::sumIndirectAndInternalFaces(mvvfOwn, mvvfNei);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        vector gradf = (Cf[facei] - C[own]) & g[own];

        scalar vsfOwn = gradf & vvfOwn[facei];
        scalar vsfNei = gradf & vvfNei[facei];
        scalar mvsfOwn = gradf & mvvfOwn[facei];
        scalar mvsfNei = gradf & mvvfNei[facei];

        scalar maxFace = max(mvsfOwn, mvsfNei);
        scalar minFace = min(mvsfOwn, mvsfNei);
        scalar maxMinFace = rk*(maxFace - minFace);
        maxFace += maxMinFace;
        minFace -= maxMinFace;

        limitFace
        (
            limiter[own],
            maxFace - vsfOwn, minFace - vsfOwn,
            magSqr(gradf)
        );


        // neighbour side
        gradf = (Cf[facei] - C[nei]) & g[nei];

        vsfOwn = gradf & vvfOwn[facei];
        vsfNei = gradf & vvfNei[facei];
        mvsfOwn = gradf & mvvfOwn[facei];
        mvsfNei = gradf & mvvfNei[facei];

        maxFace = max(mvsfOwn, mvsfNei);
        minFace = min(mvsfOwn, mvsfNei);

        limitFace
        (
            limiter[nei],
            maxFace - vsfNei, minFace - vsfNei,
            magSqr(gradf)
        );
    }


    const volVectorField::Boundary& bvf = vvf.boundaryField();

    forAll(bvf, patchi)
    {
        const fvPatchVectorField& pvf = bvf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        const vectorField& pvvfOwn(vvfOwn.boundaryField()[patchi]);
        const vectorField& pmvvfOwn(mvvfOwn.boundaryField()[patchi]);
        const vectorField& pmvvfNei(mvvfNei.boundaryField()[patchi]);

        if (pvf.coupled() || pvf.fixesValue())
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];

                vector gradf = (pCf[pFacei] - C[own]) & g[own];

                scalar vsfOwn = gradf & pvvfOwn[pFacei];
                scalar mvsfOwn = gradf & pmvvfOwn[pFacei];
                scalar mvsfNei = gradf & pmvvfNei[pFacei];

                scalar maxFace = max(mvsfOwn, mvsfNei);
                scalar minFace = min(mvsfOwn, mvsfNei);
                scalar maxMinFace = rk*(maxFace - minFace);
                maxFace += maxMinFace;
                minFace -= maxMinFace;

                limitFace
                (
                    limiter[own],
                    maxFace - vsfOwn, minFace - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vvf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.primitiveFieldRef() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// ************************************************************************* //
