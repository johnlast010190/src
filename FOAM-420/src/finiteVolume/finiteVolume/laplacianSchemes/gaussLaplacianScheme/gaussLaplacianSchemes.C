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

#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(gaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvmLaplacian         \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf              \
    (                                                                          \
        gamma*mesh.magSf()                                                     \
    );                                                                         \
                                                                               \
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected                         \
    (                                                                          \
        gammaMagSf,                                                            \
        this->tsnGradScheme_().deltaCoeffs(vf),                                \
        vf                                                                     \
    );                                                                         \
    fvMatrix<Type>& fvm = tfvm.ref();                                          \
                                                                               \
    if (this->tsnGradScheme_().corrected())                                    \
    {                                                                          \
        if (mesh.schemes().fluxRequired(vf.name()))                                      \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = new                                  \
            GeometricField<Type, fvsPatchField, surfaceMesh>                   \
            (                                                                  \
                gammaMagSf*this->tsnGradScheme_().correction(vf)               \
            );                                                                 \
                                                                               \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    *fvm.faceFluxCorrectionPtr()                               \
                )().primitiveField();                                          \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)           \
                )().primitiveField();                                          \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh>> \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian         \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian                \
    (                                                                          \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())         \
    );                                                                         \
                                                                               \
    tLaplacian.ref().rename                                                    \
    (                                                                          \
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'                    \
    );                                                                         \
                                                                               \
    return tLaplacian;                                                         \
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


#define declareFvmScalarBLaplacianGamma(Type)                                  \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvBlockMatrix<Foam::scalar>>                                   \
Foam::fv::gaussLaplacianScheme<Foam::scalar, Foam::Type>::fvmBLaplacian        \
(                                                                              \
    const GeometricField<Type, fvsPatchField, surfaceMesh>& gamma,             \
    const GeometricField<scalar, fvPatchField, volMesh>& vf                    \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());                       \
                                                                               \
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);                       \
    const GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf        \
    (                                                                          \
        SfGamma & Sn                                                           \
    );                                                                         \
                                                                               \
    tmp<surfaceScalarField> tdeltaCoeffs                                       \
    (                                                                          \
        this->tsnGradScheme_().deltaCoeffs(vf)                                 \
    );                                                                         \
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();                    \
                                                                               \
    tmp<fvBlockMatrix<scalar>> tfvbm                                           \
    (                                                                          \
        new fvBlockMatrix<scalar>                                              \
        (                                                                      \
            const_cast<GeometricField<scalar, fvPatchField, volMesh>&>(vf)     \
        )                                                                      \
    );                                                                         \
    fvBlockMatrix<scalar>& fvbm = tfvbm.ref();                                 \
    scalarField& source = fvbm.source();                                       \
                                                                               \
    CoeffField<scalar>::linearTypeField& d = fvbm.diag().asScalar();           \
    CoeffField<scalar>::linearTypeField& u = fvbm.upper().asScalar();          \
                                                                               \
    u = deltaCoeffs.primitiveField()*gammaMagSf.internalField();               \
    fvbm.negSumDiag();                                                         \
                                                                               \
    forAll(vf.boundaryField(), pI)                                             \
    {                                                                          \
        const fvPatchField<scalar>& pvf = vf.boundaryField()[pI];              \
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[pI];    \
        const fvsPatchScalarField& pDeltaCoeffs =                              \
            deltaCoeffs.boundaryField()[pI];                                   \
                                                                               \
        const fvPatchScalarField& pf = vf.boundaryField()[pI];                 \
        const fvPatch& patch = pf.patch();                                     \
        const labelUList& fc = patch.faceCells();                              \
                                                                               \
        const scalarField iCoeffs(pvf.gradientInternalCoeffs(pDeltaCoeffs));   \
                                                                               \
        forAll(pf, fI)                                                         \
        {                                                                      \
            d[fc[fI]] += pGamma[fI]*iCoeffs[fI];                               \
        }                                                                      \
                                                                               \
        if (pvf.coupled())                                                     \
        {                                                                      \
            CoeffField<scalar>::linearTypeField& pcUpper =                     \
                fvbm.coupleUpper()[pI].asScalar();                             \
                                                                               \
            pcUpper -= pGamma*iCoeffs;                                         \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            const scalarField bCoeffs(pvf.gradientBoundaryCoeffs(pDeltaCoeffs));\
                                                                               \
            forAll(pf, fI)                                                     \
            {                                                                  \
                source[fc[fI]] -= pGamma[fI]*bCoeffs[fI];                      \
            }                                                                  \
                                                                               \
        }                                                                      \
    }                                                                          \
    FatalErrorInFunction                                                       \
        << "scalar not suppoted"                                               \
        << abort(FatalError);                                                  \
                                                                               \
    return tfvbm;                                                              \
}                                                                              \


declareFvmScalarBLaplacianGamma(sphericalTensor);
declareFvmScalarBLaplacianGamma(symmTensor);
declareFvmScalarBLaplacianGamma(tensor);

// ************************************************************************* //
