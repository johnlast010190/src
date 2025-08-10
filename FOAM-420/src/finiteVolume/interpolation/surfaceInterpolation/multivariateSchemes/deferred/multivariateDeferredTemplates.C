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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/multivariateSchemes/deferred/multivariateDeferred.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchField.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::multivariateDeferred<Type>::fieldScheme::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        scheme_->corrected() ?
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "multivariateDeferred::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            scheme_->correction(vf)
        ) :
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "multivariateDeferred::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    tmp<surfaceScalarField> tWeights = scheme_->weights(vf);
    const surfaceScalarField& weights = tWeights();

    // Correction is scheme value minus upwind (scheme correction added above)
    forAll(faceFlux, facei)
    {
        scalar fw(weights[facei]);
        label celli = (faceFlux[facei] >= 0) ? owner[facei] : neighbour[facei];
        sfCorr[facei]
            = fw*vf[owner[facei]] + (1-fw)*vf[neighbour[facei]] - vf[celli];
    }


    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary&
        bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner =
                mesh.boundary()[patchi].faceCells();

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            tmp<Field<Type>> tvfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<Type>& vfNei = tvfNei();

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                scalar fw(weights.boundaryField()[patchi][facei]);
                pSfCorr[facei] = fw*vf[own] + (1-fw)*vfNei[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] -= vf[own];
                }
                else
                {
                    pSfCorr[facei] -= vfNei[facei];
                }
            }
        }
        else
        {
            pSfCorr = pTraits<Type>::zero;
        }
    }

    if (scheme_->corrected())
    {
        sfCorr += scheme_->correction(vf);
    }

    return tsfCorr;
}


// ************************************************************************* //
