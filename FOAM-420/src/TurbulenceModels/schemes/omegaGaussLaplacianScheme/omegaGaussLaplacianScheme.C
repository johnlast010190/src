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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "omegaGaussLaplacianScheme.H"
#include "fvMesh/fvMesh.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "turbulenceModel.H"
#include "derivedFvPatchFields/wallFunctions/omegaWallFunctions/knopOmegaWallFunction/knopOmegaWallFunctionFvPatchScalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Generalise dot product to avoid unnecessary specialisations
tmp<surfaceVectorField> operator&
(
    const surfaceVectorField& svf, const surfaceScalarField& ssf
)
{
    return svf * ssf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace fv
{

template<class Type, class GType>
scalar omegaGaussLaplacianScheme<Type, GType>::getDeltaStarCoefficient
(
    const label& patchID,
    const label& patchFaceLabel,
    const label& faceLabel,
    const scalar& yPlus
)
{
    scalar deltaStar = 1;
    const fvMesh& mesh = this->mesh();

    const fvPatch& patch = mesh.boundary()[patchID];
    const vectorField& Sf(patch.Sf());
    const scalarField& magSf(patch.magSf());
    const vectorField& Cf(patch.Cf());

    const label& faceOwnerLabel = mesh.owner()[faceLabel];
    const label& faceNeihborLabel = mesh.neighbour()[faceLabel];

    const point& cellCentreOwn = mesh.C()[faceOwnerLabel];
    const point& cellCentreNei = mesh.C()[faceNeihborLabel];
    const vector& boundaryFaceNormal = Sf[patchFaceLabel]/magSf[patchFaceLabel];
    const vector& currentFaceNormal =
        mesh.Sf()[faceLabel] / (mesh.magSf()[faceLabel] + VSMALL);
    scalar weight =
        min
            (
                1.,
                mag(boundaryFaceNormal & currentFaceNormal)
            );

    scalar distOwn = mag((cellCentreOwn - Cf[patchFaceLabel])&currentFaceNormal);
    scalar distNei = mag((cellCentreNei - Cf[patchFaceLabel])&currentFaceNormal);
    scalar distFace = mag((mesh.Cf()[faceLabel] - Cf[patchFaceLabel])&currentFaceNormal);

    //Limit Minimum Accepted Distances
    distOwn = max(ROOTSMALL, distOwn);
    distNei = max(ROOTSMALL, distNei);
    distFace = max(ROOTSMALL, distFace);

    if (yPlus<3.)
    {
        scalar w1 = 1./(1.+exp(2.*(yPlus-2.)));
        scalar w2 = 1./(1.+exp(-2.*(yPlus-2.)));
        scalar deltaStar_upper = sqr(2.*sqr(distNei)*sqr(distOwn)/
                ((distOwn+distNei)*pow3(distFace)));
        scalar deltaStar_lower = 2.*sqr(distNei)*sqr(distOwn)/
                ((distOwn+distNei)*pow3(distFace));
        deltaStar = deltaStar_upper*w1 + deltaStar_lower*w2;
    }

    if (yPlus >= 3. && yPlus < 10.)
    {
        deltaStar = 2.*sqr(distNei)*sqr(distOwn)/
            ((distOwn+distNei)*pow3(distFace));

//        deltaStar = distOwn*distNei/sqr(distFace);
        scalar w =  1/cosh(1.8*(yPlus-3.));
        deltaStar = pow(deltaStar, w);
    }
    deltaStar = pow(deltaStar, weight);
    deltaStar = min(max(0.3, deltaStar), 1.);

    if (debug)
    {
        Info<< "y+ is " << yPlus << endl;
        Info<< "Delta Star " << deltaStar << endl;
        Info<< "Weight of face " << faceLabel << ": " << weight
             << endl;
    }
    return deltaStar;
}

template<class Type, class GType>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fv::omegaGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions() * gammaMagSf.dimensions()
                * vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    tmp<surfaceScalarField> tmGammaMagSf = fvc::applyFaceMask(gammaMagSf);
    const surfaceScalarField& mGammaMagSf = tmGammaMagSf();


    fvm.upper() = deltaCoeffs.primitiveField() * mGammaMagSf.internalField();

    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const fvsPatchScalarField
            & pmGamma = mGammaMagSf.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];

        if (pvf.coupled())
        {
            fvm.internalCoeffs()[patchi] =
                pmGamma * pvf.gradientInternalCoeffs(pDeltaCoeffs);
            fvm.boundaryCoeffs()[patchi] =
                -pmGamma * pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
        }
        else
        {
            fvm.internalCoeffs()[patchi] =
                pmGamma * pvf.gradientInternalCoeffs();
            fvm.boundaryCoeffs()[patchi] =
                -pmGamma * pvf.gradientBoundaryCoeffs();
        }
    }

    return tfvm;
}

template<class Type, class GType>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::omegaGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceVectorField& SfGammaCorr,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                        (
                            "gammaSnGradCorr(" + vf.name() + ')',
                            vf.instance(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                    mesh,
                    SfGammaCorr.dimensions()
                        * vf.dimensions() * mesh.deltaCoeffs().dimensions()
                )
        );
    tgammaSnGradCorr.ref().oriented() = SfGammaCorr.oriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tgammaSnGradCorr.ref().replace
            (
                cmpt,
                fvc::dotInterpolate(SfGammaCorr, fvc::grad(vf.component(cmpt)))
            );
    }

    return tgammaSnGradCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fv::omegaGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
        (
            fvc::div(this->tsnGradScheme_().snGrad(vf) * mesh.magSf())
        );

    tLaplacian.ref().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}

template<class Type, class GType>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fv::omegaGaussLaplacianScheme<Type, GType>::fvmLaplacian
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    )
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf() / mesh.magSf());

    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
        (
            SfGamma & Sn
        );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn * Sn);

    surfaceScalarField deltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const turbulenceModel& turb = vf.db().template lookupObject<turbulenceModel>
        (
            IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    vf.internalField().group()
                )
        );

    List<bool> isVisited(mesh.faces().size(), false);
    const auto& omega = vf.db().template lookupObject<volScalarField>("omega");
    const volScalarField::Boundary& bf = omega.boundaryField();
    forAll(bf, patchi)
    {
        if
        (
            isA<knopOmegaWallFunctionFvPatchScalarField>(bf[patchi])
        )
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            const tmp<volScalarField> tnut = turb.nut();
            const volScalarField& nut = tnut();
            const fvPatchScalarField& nutwpf = nut.boundaryField()[patchi];
            const auto& nutWF =
                dynamic_cast<const nutWallFunctionFvPatchScalarField&>(nutwpf);
            const scalarField& omegaWall = omega.boundaryField()[patchi];
            tmp<scalarField> yp  = nutWF.yPlus();
            forAll(omegaWall, fi)
            {
                const scalar& yPlus = yp()[fi];
                label faceCellI = patch.faceCells()[fi];
                const cell& fCell = mesh.cells()[faceCellI];
                forAll(fCell, ffi)
                {
                    const label& facei = fCell[ffi];
                    if (facei < mesh.nInternalFaces() && !isVisited[facei])
                    {
                        isVisited[facei] = true;
                        scalar deltaStar = getDeltaStarCoefficient
                        (
                            patchi, fi, facei, yPlus
                        );

                        deltaCoeffs[facei] *= deltaStar;
                    }
                }
            }
        }
    }
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected
    (
        SfGammaSn,
        deltaCoeffs,
        vf
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfaceFluxCorrection
        = gammaSnGradCorr(SfGammaCorr, vf);

    if (this->tsnGradScheme_().corrected())
    {
        tfaceFluxCorrection.ref() +=
            SfGammaSn * this->tsnGradScheme_().correction(vf);
    }

    fvm.source() -=
        mesh.V() * fvc::div(tfaceFluxCorrection())().primitiveField();

    if (mesh.schemes().fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
    }

    return tfvm;
}

template<class Type, class GType>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fv::omegaGaussLaplacianScheme<Type, GType>::fvcLaplacian
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    )
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf() / mesh.magSf());
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
        (
            SfGamma & Sn
        );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn * Sn);

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
        (
            fvc::div
                (
                    SfGammaSn * this->tsnGradScheme_().snGrad(vf)
                        + gammaSnGradCorr(SfGammaCorr, vf)
                )
        );

    tLaplacian.ref().rename
        (
            "laplacian(" + gamma.name() + ',' + vf.name() + ')'
        );

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
// ************************************************************************* //
