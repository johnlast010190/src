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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchFields.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc, class BlendingFactor>
void Foam::LimitedICScheme<Type, Limiter, LimitFunc, BlendingFactor>::calcLimiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    surfaceScalarField& limiterField
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh>>
        tlPhi = LimitFunc<Type>()(phi);

    const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
        lPhi = tlPhi();

    tmp<GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>>
        tgradc(fvc::grad(lPhi));
    const GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>&
        gradc = tgradc();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = limiterField.primitiveFieldRef();

    // calculate the face Courant number
    const surfaceScalarField faceFlux = this->faceFlux_;
    volScalarField outFlux ( fvc::surfaceSum(mag(faceFlux)) );
    outFlux.primitiveFieldRef() /= mesh.V(); //internalField

    //- fix mass conservation problems with coupled boundaries (parallel runs)
    //  i.e. correct calculation of local Courant number at coupled boundaries
    if (Pstream::parRun())
    {
        outFlux.correctBoundaryConditions();
    }

    surfaceScalarField Cof
    (
        0.5 * mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux).interpolate
        (
            outFlux
        )
    );

    // get indicator field from objectRegistry based on the fieldName
    const volScalarField& vsf =
    mesh.objectRegistry::lookupObject<const volScalarField>
        (
            fieldName_
        );

    // calculate the blending factors
    tmp<surfaceScalarField> tblendingFactor = BlendingFactor()(vsf);
    surfaceScalarField& blendingFactor = tblendingFactor.ref();

    // for testing: write out the blending field
    //volScalarField blending = fvc::average(blendingFactor);
    //surfaceScalarField blendingf = blendingFactor;
    //blending.write();
    //blendingf.write();
    blendingFactor = blendingFactor*0. +1;

    forAll(pLim, face)
    {
        label own = owner[face];
        label nei = neighbour[face];

        point p;

        if (this->faceFlux_[face] > 0)
        {
            p = C[own] - (C[nei] - C[own]);
        }
        else
        {
            p = C[nei] - (C[own] - C[nei]);
        }

        pLim[face] = Limiter::limiter
        (
            CDweights[face],
            this->faceFlux_[face],
            lPhi[own],
            lPhi[nei],
            gradc[own],
            gradc[nei],
            C[nei] - C[own],
            Cof[face],
            phi,
            p,
            face,
            blendingFactor[face]
        );
    }

    surfaceScalarField::Boundary& bLim =
        limiterField.boundaryFieldRef();

    forAll(bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            const Field<typename Limiter::phiType> plPhiP
            (
                lPhi.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::phiType> plPhiN
            (
                lPhi.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<typename Limiter::gradPhiType> pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::gradPhiType> pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(CDweights.boundaryField()[patchi].patch().delta());

            // Get the Courant number on the boundary patches
            const scalarField& pCof = Cof.boundaryField()[patchi];

            // Get blending factors on the boundary patches
            const scalarField& pblendingFactor = blendingFactor.boundaryField()[patchi];

            forAll(pLim, face)
            {
                pLim[face] = Limiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    plPhiP[face],
                    plPhiN[face],
                    pGradcP[face],
                    pGradcN[face],
                    pd[face],
                    pCof[face],
                    phi,
                    pTraits<vector>::zero,
                    face,
                    pblendingFactor[face]
                );
            }
        }
        else
        {
            pLim = 1.0;
        }
    }
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc, class BlendingFactor>
Foam::tmp<Foam::surfaceScalarField>
Foam::LimitedICScheme<Type, Limiter, LimitFunc, BlendingFactor>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    const word limiterFieldName(type() + "Limiter(" + phi.name() + ')');

    if (this->mesh().solution().cache("limiter"))
    {
        if (!mesh.foundObject<surfaceScalarField>(limiterFieldName))
        {
            surfaceScalarField* limiterField
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        limiterFieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimless
                )
            );

            mesh.objectRegistry::store(limiterField);
        }

        surfaceScalarField& limiterField =
            const_cast<surfaceScalarField&>
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    limiterFieldName
                )
            );

        calcLimiter(phi, limiterField);

        return limiterField;
    }
    else
    {
        tmp<surfaceScalarField> tlimiterField
        (
            new surfaceScalarField
            (
                IOobject
                (
                    limiterFieldName,
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimless
            )
        );

        calcLimiter(phi, tlimiterField.ref());

        return tlimiterField;
    }
}

// ************************************************************************* //
