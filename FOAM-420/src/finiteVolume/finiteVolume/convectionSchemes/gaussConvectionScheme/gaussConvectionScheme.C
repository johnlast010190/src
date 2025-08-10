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

#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const surfaceInterpolationScheme<Type>&
gaussConvectionScheme<Type>::interpScheme() const
{
    return tinterpScheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussConvectionScheme<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
gaussConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<surfaceScalarField> tmFaceFlux = fvc::applyFaceMask(faceFlux);
    const surfaceScalarField& mFaceFlux = tmFaceFlux();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.lower() = -weights.primitiveField()*mFaceFlux.primitiveField();
    fvm.upper() = fvm.lower() + mFaceFlux.primitiveField();

    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& mPatchFlux =
            mFaceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = mPatchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -mPatchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tmFluxCorr =
            mFaceFlux*tinterpScheme_().correction(vf);
        fvm += fvc::surfaceIntegrate(tmFluxCorr());
        if (vf.mesh().schemes().fluxRequired(vf.name()))
        {
            fvm.faceFluxCorrectionPtr() = tmFluxCorr.ptr();
        }
    }

    return tfvm;
}


template<class Type>
tmp<fvBlockMatrix<Type>>
gaussConvectionScheme<Type>::fvmBDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    tmp<fvBlockMatrix<Type>> tbfvm
    (
        new fvBlockMatrix<Type>
        (
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>(vf)
        )
    );

    tmp<fvMatrix<Type>> segConv = fvmDiv(faceFlux, vf);

    tbfvm->insertEquation(0, segConv.ref());

    return tbfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
gaussConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const volScalarField& vfMult,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<surfaceScalarField> tmFaceFlux = fvc::applyFaceMask(faceFlux);
    const surfaceScalarField& mFaceFlux = tmFaceFlux();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vfMult.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.lower() = -weights.primitiveField()*mFaceFlux.primitiveField();
    fvm.upper() = fvm.lower() + mFaceFlux.primitiveField();

    const labelUList& owner = vf.mesh().owner();
    const labelUList& neighbour = vf.mesh().neighbour();
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        fvm.lower()[facei] *= vfMult[own];
        fvm.upper()[facei] *= vfMult[nei];
    }

    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& mPatchFlux =
            mFaceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] =
            mPatchFlux*vfMult.boundaryField()[patchi].patchInternalField()
           *psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] =
           -mPatchFlux*vfMult.boundaryField()[patchi]
           *psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tmFluxCorr =
            mFaceFlux*tinterpScheme_().correction(vfMult*vf);
        fvm += fvc::surfaceIntegrate(tmFluxCorr());
        if (vf.mesh().schemes().fluxRequired(vf.name()))
        {
            fvm.faceFluxCorrectionPtr() = tmFluxCorr.ptr();
        }
    }

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
gaussConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
