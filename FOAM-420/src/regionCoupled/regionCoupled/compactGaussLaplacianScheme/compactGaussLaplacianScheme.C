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
    (c) 2019-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "compactGaussLaplacianScheme/compactGaussLaplacianScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/snGradSchemes/orthogonalSnGrad/orthogonalSnGrad.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "finiteVolume/snGradSchemes/limitedSnGrad/limitedSnGrad.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generalise dot product to avoid unnecessary specialisations
tmp<surfaceVectorField> operator&
(
    surfaceVectorField svf, surfaceScalarField ssf
)
{
    return svf * ssf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

template<class Type>
tmp<surfaceScalarField> overRelaxedDeltaCoeffs
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return gamma.mesh().nonOrthDeltaCoeffs();
}

template<class Type, class TensorType>
tmp<surfaceScalarField> overRelaxedDeltaCoeffs
(
    const GeometricField<TensorType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    // Re-do the calculation for mesh nonOrthDeltaCoeffs with effective face
    // normal Sf & gamma
    const fvMesh& mesh = gamma.mesh();
    surfaceVectorField unitSfEff(mesh.Sf() & gamma);
    dimensionedScalar sm("small", unitSfEff.dimensions(), VSMALL);
    unitSfEff /= stabilise(mag(unitSfEff), sm);

    tmp<surfaceVectorField> tdelta = mesh.delta();
    tmp<surfaceScalarField> tCoeffs =
        scalar(1)/max(unitSfEff & tdelta(), 0.05*mag(tdelta()));
    tCoeffs->setOriented();

    return tCoeffs;
}


// * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
compactGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
    const surfaceVectorField lHat( mesh.delta()/mag(mesh.delta()) );

    const surfaceVectorField SfGamma(mesh.Sf());
    const surfaceScalarField SfGammaNorm( SfGamma & lHat );
    const surfaceVectorField SfGammaTang( SfGamma - SfGammaNorm*lHat );

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            mesh.deltaCoeffs().dimensions()
           *SfGammaNorm.dimensions()
           *vf.dimensions()
           /dimVolume,
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsg =
        SfGammaNorm*this->tsnGradScheme_().snGrad(vf);
    for (direction cmp = 0; cmp < pTraits<Type>::nComponents; cmp++)
    {
        tLaplacian->replace
        (
            cmp,
            fvc::div
            (
                tsg->component(cmp)
              + (SfGammaTang & linearInterpolate(fvc::grad(vf.component(cmp))))
            )()
        );
    }

    return tLaplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
compactGaussLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

    surfaceVectorField delta( mesh.delta() );
    const surfaceVectorField lHat( delta/mag(delta) );

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    const surfaceVectorField SfGamma(mSf & gamma);
    // This uses over-relaxation proportional to effective face non-orthogonality
    const surfaceScalarField SfGammaCoeff
    (
        mag(SfGamma)*mag(delta)*overRelaxedDeltaCoeffs(gamma, vf)
    );
    const surfaceScalarField SfGammaNorm( SfGamma & lHat );
    const surfaceVectorField SfGammaTang( SfGamma - SfGammaNorm*lHat );

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs().dimensions()*SfGamma.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tFluxCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            "fluxCorr",
            vf.db(),
            mesh,
            fvm.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& fluxCorr =
        tFluxCorr.ref();
    fluxCorr.setOriented();

    forAll(fvm.upper(), facei)
    {
        label own = mesh.owner()[facei];
        label nei = mesh.neighbour()[facei];

        fvm.upper()[facei] = deltaCoeffs[facei]*SfGammaCoeff[facei];
        fvm.diag()[own] += -deltaCoeffs[facei]*SfGammaCoeff[facei];

        // Keep matrix symmetric
        //fvm.lower()[facei] =
        //    deltaCoeffs[facei]*SfGammaCoeff[facei];
        fvm.diag()[nei] += -deltaCoeffs[facei]*SfGammaCoeff[facei];

        // Remove the above; replace with correction
        fluxCorr[facei] =
            deltaCoeffs[facei]
           *(SfGammaNorm[facei]-SfGammaCoeff[facei])
           *(vf[nei]-vf[own]);
    }

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& pSfGammaCoeff =
            SfGammaCoeff.boundaryField()[patchi];
        const fvsPatchScalarField& pSfGammaNorm =
            SfGammaNorm.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];
        fvsPatchField<Type>& pFluxCorr = fluxCorr.boundaryFieldRef()[patchi];

        const Field<Type> pivf( pvf.patchInternalField() );
        if (pvf.coupled())
        {
            const Field<Type> pnvf( pvf.patchNeighbourField() );

            fvm.boundaryCoeffs()[patchi] =
               -pSfGammaCoeff*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            fvm.internalCoeffs()[patchi] =
                pSfGammaCoeff*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            // Remove the above; replace with correction
            for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                const scalarField ps
                (
                    (pSfGammaNorm-pSfGammaCoeff)
                   *pvf.gradientBoundaryCoeffs(pDeltaCoeffs)().component(cmpt)
                   *pnvf.component(cmpt)
                  + (pSfGammaNorm-pSfGammaCoeff)
                   *pvf.gradientInternalCoeffs(pDeltaCoeffs)().component(cmpt)
                   *pivf.component(cmpt)
                );
                pFluxCorr.replace(cmpt, ps);
            }
        }
        else
        {
            fvm.boundaryCoeffs()[patchi] =
               -pSfGammaCoeff*pvf.gradientBoundaryCoeffs();
            fvm.internalCoeffs()[patchi] =
                pSfGammaCoeff*pvf.gradientInternalCoeffs();
            if (pvf.regionCoupled())
            {
                const regionCoupledFvPatchField<Type>& rcfvpf =
                    refCast<const regionCoupledFvPatchField<Type>>(pvf);
                // Remove the above; replace with correction
                for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
                {
                    const scalarField ps
                    (
                        (pSfGammaNorm-pSfGammaCoeff)
                       *pvf.gradientBoundaryCoeffs()().component(cmpt)
                       *(
                           pvf.component(cmpt)
                           // Don't over-relax the explicit correction
                         - rcfvpf.faceCorrEval()().component(cmpt)
                         + rcfvpf.faceCorr()().component(cmpt)
                        )
                      + (pSfGammaNorm-pSfGammaCoeff)
                       *pvf.gradientInternalCoeffs()().component(cmpt)
                       *pivf.component(cmpt)
                    );
                    pFluxCorr.replace(cmpt, ps);
                }
            }
            else
            {
                // Remove the above; replace with correction
                for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
                {
                    const scalarField ps
                    (
                        (pSfGammaNorm-pSfGammaCoeff)
                       *pvf.gradientBoundaryCoeffs()().component(cmpt)
                      + (pSfGammaNorm-pSfGammaCoeff)
                       *pvf.gradientInternalCoeffs()().component(cmpt)
                       *pivf.component(cmpt)
                    );
                    pFluxCorr.replace(cmpt, ps);
                }
            }
        }
    }

    autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> tangCorr;
    if (isA<limitedSnGrad<Type>>(this->tsnGradScheme_()))
    {
        tangCorr.set
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                "tangCorr",
                vf.db(),
                mesh,
                vf.dimensions()/dimLength
            )
        );
    }

    for (direction cmp = 0; cmp < pTraits<Type>::nComponents; cmp++)
    {
        // If possible, use precalculated grad (even if not cached - for sharing
        // with BCs)
        tmp<volScalarField> tvfCmp = vf.component(cmp);
        word gradName = "grad(" + tvfCmp().name() + ")";

        surfaceVectorField grad
        (
            linearInterpolate
            (
                mesh.foundObject<volVectorField>(gradName) ?
                mesh.lookupObject<volVectorField>(gradName) :
                fvc::grad(tvfCmp)()
            )
        );

        fluxCorr.replace
        (
            cmp,
            fluxCorr.component(cmp) + (SfGammaTang & grad)
        );

        if (tangCorr.valid())
        {
            const surfaceVectorField SfHat( mSf/mesh.magSf() );
            const surfaceVectorField SfHatTang( SfHat - (SfHat&lHat)*lHat );
            tangCorr->replace(cmp, (SfHatTang & grad));
        }
    }

    if (isA<limitedSnGrad<Type>>(this->tsnGradScheme_()))
    {
        scalar limitCoeff =
            refCast<const limitedSnGrad<Type>>(this->tsnGradScheme_()).limitCoeff();

        // Calcuate limiting as if there was no gamma - since the effect of
        // uniform non-orthogonality of the effective faces (Sf&gamma) cancels
        // when taking the div of the explicit correction, and we assume gamma
        // does not change rapidly

        GeometricField<Type, fvsPatchField, surfaceMesh> snGradImpl
        (
            snGradScheme<Type>::snGrad(vf, mesh.nonOrthDeltaCoeffs(), "SndGrad")
        );
        GeometricField<Type, fvsPatchField, surfaceMesh> snGradExpl
        (
            snGradScheme<Type>::snGrad(vf, mesh.deltaCoeffs(), "SndGrad")
        );
        GeometricField<Type, fvsPatchField, surfaceMesh> corr
        (
            snGradExpl - snGradImpl + tangCorr()
        );

        surfaceScalarField limiter
        (
            min
            (
                limitCoeff
               *mag(snGradImpl)
               /(
                    (1 - limitCoeff)*mag(corr)
                  + dimensionedScalar("small", corr.dimensions(), SMALL)
                ),
                dimensionedScalar("one", dimless, 1.0)
            )
        );

        // Don't limit ordinary boundaries as the gradient can be constrained
        // by the BC and this leads to incorrect limiting of the
        // explicit part in anisotropic cases
        forAll(limiter.boundaryField(), patchi)
        {
            if
            (
                !vf.boundaryField()[patchi].coupled()
             && !vf.boundaryField()[patchi].regionCoupled()
            )
            {
                limiter.boundaryFieldRef()[0] = 1.0;
            }
        }

        fluxCorr *= limiter;
    }

    // Goes on RHS
    fvm.source() = -fvc::div(fluxCorr)().primitiveField()*mesh.Vsc();

    if (mesh.schemes().fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tFluxCorr.ptr();
    }

    return tfvm;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
compactGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
    const surfaceVectorField lHat( mesh.delta()/mag(mesh.delta()) );

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    const surfaceVectorField SfGamma(mSf & gamma);
    const surfaceScalarField SfGammaNorm( SfGamma & lHat );
    const surfaceVectorField SfGammaTang( SfGamma - SfGammaNorm*lHat );

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            mesh.deltaCoeffs().dimensions()
           *SfGammaNorm.dimensions()
           *vf.dimensions()
           /dimVolume,
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
    for (direction cmp = 0; cmp < pTraits<Type>::nComponents; cmp++)
    {
        tLaplacian->replace
        (
            cmp,
            fvc::div
            (
                SfGammaNorm
               *fv::orthogonalSnGrad<scalar>(mesh).snGrad(vf.component(cmp))
              + (SfGammaTang & linearInterpolate(fvc::grad(vf.component(cmp))))
            )()
        );
    }

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
