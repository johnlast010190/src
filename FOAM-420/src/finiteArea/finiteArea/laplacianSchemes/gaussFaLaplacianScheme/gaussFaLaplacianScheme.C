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

#include "finiteArea/laplacianSchemes/gaussFaLaplacianScheme/gaussFaLaplacianScheme.H"
#include "interpolation/edgeInterpolation/edgeInterpolate.H"
#include "finiteArea/fac/facDiv.H"
#include "finiteArea/fac/facGrad.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type, class GType>
tmp<faMatrix<Type>>
gaussLaplacianScheme<Type, GType>::famLaplacianUncorrected
(
    const edgeScalarField& gammaMagSf,
    const edgeScalarField& deltaCoeffs,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
    fam.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const faPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const faePatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
        const faePatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];

        if (pvf.coupled())
        {
            fam.internalCoeffs()[patchi] =
                pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            fam.boundaryCoeffs()[patchi] =
               -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
        }
        else
        {
            fam.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
            fam.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
        }
    }

    return tfam;
}

template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
gaussLaplacianScheme<Type, GType>::facLaplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = this->mesh();

    const edgeVectorField Sn(mesh.edgeAreaNormals());
    const edgeVectorField SfGamma(mesh.Le() & gamma);
    const GeometricField<scalar, faePatchField, faEdgeMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const edgeVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

    tmp<GeometricField<Type, faPatchField, areaMesh>> tLaplacian
    (
        fac::div
        (
            SfGammaSn*this->tlnGradScheme_().lnGrad(vf)
          + gammaLnGradCorr(SfGammaCorr, vf)
        )
    );

    tLaplacian.ref().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

    return tLaplacian;
}

template<class Type, class GType>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
gaussLaplacianScheme<Type, GType>::gammaLnGradCorr
(
    const edgeVectorField& SfGammaCorr,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = this->mesh();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tgammaLnGradCorr
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "gammaLnGradCorr("+vf.name()+')',
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tgammaLnGradCorr.ref().replace
        (
            cmpt,
            SfGammaCorr & fac::interpolate(fac::grad(vf.component(cmpt)))
//            fac::interpolate(fac::grad(vf.component(cmpt)))
        );
    }

    return tgammaLnGradCorr;
}

template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
gaussLaplacianScheme<Type, GType>::facLaplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = this->mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh>> tLaplacian
    (
        fac::div(this->tlnGradScheme_().lnGrad(vf)*mesh.magLe())
    );

    tLaplacian.ref().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}

template<class Type, class GType>
tmp<faMatrix<Type>>
gaussLaplacianScheme<Type, GType>::famLaplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = this->mesh();

    const edgeVectorField En(mesh.Le()/mesh.magLe());

    const edgeVectorField EfGamma(mesh.Le() & gamma);
    const GeometricField<scalar, faePatchField, faEdgeMesh> EfGammaEn
    (
        EfGamma & En
    );
    const edgeVectorField EfGammaCorr(EfGamma - EfGammaEn*En);

    tmp<faMatrix<Type>> tfam = famLaplacianUncorrected
    (
        EfGammaEn,
        this->tlnGradScheme_().deltaCoeffs(vf),
        vf
    );
    faMatrix<Type>& fam = tfam.ref();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tfaceFluxCorrection
        = gammaLnGradCorr(EfGammaCorr, vf);

    if (this->tlnGradScheme_().corrected())
    {
        tfaceFluxCorrection.ref() +=
            EfGammaEn*this->tlnGradScheme_().correction(vf);
    }

    fam.source() -= mesh.S()*fac::div(tfaceFluxCorrection())().internalField();

    if (mesh.schemes().fluxRequired(vf.name()))
    {
        fam.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
    }

    return tfam;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
