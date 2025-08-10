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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/schemes/linearUpwind/linearUpwindV.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace
{

template <class Type>
Foam::scalar luInnerProduct (Type v1, Type v2)
{
    return v1 & v2;
}


template <>
Foam::scalar luInnerProduct (Foam::scalar s1, Foam::scalar s2)
{
    return s1*s2;
}


}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::linearUpwindV<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "linearUpwindV::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>
            (
                vf.name(),
                vf.dimensions(),
                Zero
            ),
            fvsPatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();
    const surfaceScalarField& faceFlux = this->faceFlux_;
    const surfaceScalarField& w = mesh.weights();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > tgradVf = gradScheme_().grad(vf, gradSchemeName_);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gradVf = tgradVf();

    forAll(faceFlux, facei)
    {
        Type maxCorr;

        if (faceFlux[facei] > 0.0)
        {
            maxCorr =
                (1.0 - w[facei])*(vf[nei[facei]] - vf[own[facei]]);

            sfCorr[facei] =
                (Cf[facei] - C[own[facei]]) & gradVf[own[facei]];
        }
        else
        {
            maxCorr =
                w[facei]*(vf[own[facei]] - vf[nei[facei]]);

            sfCorr[facei] =
               (Cf[facei] - C[nei[facei]]) & gradVf[nei[facei]];
        }

        scalar sfCorrs = magSqr(sfCorr[facei]);
        scalar maxCorrs = luInnerProduct(sfCorr[facei], maxCorr);

        if (maxCorrs < 0)
        {
            sfCorr[facei] = Zero;
        }
        else if (sfCorrs > maxCorrs)
        {
            sfCorr[facei] *= maxCorrs/(sfCorrs + VSMALL);
        }
        /*
        if (sfCorrs > 0)
        {
            if (maxCorrs < 0)
            {
                sfCorr[facei] = Zero;
            }
            else if (sfCorrs > maxCorrs)
            {
                sfCorr[facei] *= maxCorrs/(sfCorrs + VSMALL);
            }
        }
        else if (sfCorrs < 0)
        {
            if (maxCorrs > 0)
            {
                sfCorr[facei] = Zero;
            }
            else if (sfCorrs < maxCorrs)
            {
                sfCorr[facei] *= maxCorrs/(sfCorrs - VSMALL);
            }
        }
        */
    }


    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner =
                mesh.boundary()[patchi].faceCells();

            const vectorField& pCf = Cf.boundaryField()[patchi];
            const scalarField& pW = w.boundaryField()[patchi];

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const Field<typename outerProduct<vector, Type>::type> pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            const Field<Type> pVfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                Type maxCorr;

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];

                    maxCorr = (1.0 - pW[facei])*(pVfNei[facei] - vf[own]);
                }
                else
                {
                    pSfCorr[facei] =
                        (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];

                    maxCorr = pW[facei]*(vf[own] - pVfNei[facei]);
                }

                scalar pSfCorrs = magSqr(pSfCorr[facei]);
                scalar maxCorrs = luInnerProduct(pSfCorr[facei], maxCorr);

                if (pSfCorrs > 0)
                {
                    if (maxCorrs < 0)
                    {
                        pSfCorr[facei] = Zero;
                    }
                    else if (pSfCorrs > maxCorrs)
                    {
                        pSfCorr[facei] *= maxCorrs/(pSfCorrs + VSMALL);
                    }
                }
            }
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makelimitedSurfaceInterpolationTypeScheme(linearUpwindV, scalar)
    makelimitedSurfaceInterpolationTypeScheme(linearUpwindV, vector)
}


// ************************************************************************* //
