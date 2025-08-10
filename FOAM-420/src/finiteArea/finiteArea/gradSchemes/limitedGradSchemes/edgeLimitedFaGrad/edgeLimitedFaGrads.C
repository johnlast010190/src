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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteArea/gradSchemes/limitedGradSchemes/edgeLimitedFaGrad/edgeLimitedFaGrad.H"
#include "finiteArea/gradSchemes/gaussFaGrad/gaussFaGrad.H"
#include "faMesh/faMesh.H"
#include "areaMesh/areaMesh.H"
#include "faEdgeMesh/faEdgeMesh.H"
#include "fields/areaFields/areaFields.H"
#include "fields/faPatchFields/basic/fixedValue/fixedValueFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
makeFaGradScheme(edgeLimitedGrad)

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<areaVectorField> edgeLimitedGrad<scalar>::grad
(
    const areaScalarField& vsf
) const
{
    const faMesh& mesh = vsf.mesh();

    tmp<areaVectorField> tGrad = basicGradScheme_().grad(vsf);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaVectorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        label own = owner[edgei];
        label nei = neighbour[edgei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        // owner side
        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn, minEdge - vsfOwn,
            (Cf[edgei] - C[own]) & g[own]
        );

        // neighbour side
        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei, minEdge - vsfNei,
            (Cf[edgei] - C[nei]) & g[nei]
        );
    }

    const areaScalarField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const faPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            scalarField psfNei ( psf.patchNeighbourField() );

            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psfNei[pEdgei];

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    (pCf[pEdgei] - C[own]) & g[own]
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psf[pEdgei];

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    (pCf[pEdgei] - C[own]) & g[own]
                );
            }
        }
    }

    if (fa::debug)
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
tmp<areaTensorField> edgeLimitedGrad<vector>::grad
(
    const areaVectorField& vvf
) const
{
    const faMesh& mesh = vvf.mesh();

    tmp<areaTensorField> tGrad = basicGradScheme_().grad(vvf);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // create limiter
    scalarField limiter(vvf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        label own = owner[edgei];
        label nei = neighbour[edgei];

        vector vvfOwn = vvf[own];
        vector vvfNei = vvf[nei];

        // owner side
        vector gradf = (Cf[edgei] - C[own]) & g[own];

        scalar vsfOwn = gradf & vvfOwn;
        scalar vsfNei = gradf & vvfNei;

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn, minEdge - vsfOwn,
            magSqr(gradf)
        );


        // neighbour side
        gradf = (Cf[edgei] - C[nei]) & g[nei];

        vsfOwn = gradf & vvfOwn;
        vsfNei = gradf & vvfNei;

        maxEdge = max(vsfOwn, vsfNei);
        minEdge = min(vsfOwn, vsfNei);

        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei, minEdge - vsfNei,
            magSqr(gradf)
        );
    }
    const areaVectorField::Boundary& bvf = vvf.boundaryField();

    forAll(bvf, patchi)
    {
        const faPatchVectorField& psf = bvf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            vectorField psfNei ( psf.patchNeighbourField() );

            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psfNei[pEdgei];

                vector gradf = (pCf[pEdgei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psf[pEdgei];

                vector gradf = (pCf[pEdgei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
    }

    if (fa::debug)
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
