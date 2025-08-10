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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "alphaContactAngle/alphaContactAngle/alphaContactAngleFvPatchScalarField.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvc/fvcSnGrad.H"
#include "finiteVolume/fvc/fvcAverage.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

// To enable the use of debug switch
namespace Foam {
    defineTypeNameAndDebug(interfaceProperties, 0);
}

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    //autoPtr<volVectorField> gradAlphaPtr;
    if (nSmoothAlpha() > 0)
    {
        Info<<"Smoothing alpha"<<endl;
        volScalarField aveAlpha =  alpha1_;
        for (int nSmooth = 0; nSmooth < nSmoothAlpha(); nSmooth++)
        {
            aveAlpha = min(aveAlpha,0.95*aveAlpha + 0.05*fvc::average(aveAlpha));
        }
        gradAlphaPtr_.reset(new volVectorField(fvc::grad(aveAlpha, "nHat")));
    }
    else
    {
        gradAlphaPtr_.reset(new volVectorField(fvc::grad(alpha1_, "nHat")));
    }

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlphaPtr_()));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


void Foam::interfaceProperties::calculateCAlpha()
{
    if (cAlphaModel_ == "dynamic")
    {
        if (debug)
        {
            Info<< "Calculating dynamic compression coefficient cAlpha. " << endl;
        }

        const fvMesh& mesh = alpha1_.mesh();
        const vectorField& C = mesh.C();
        //const surfaceVectorField& Sf = mesh.Sf();
        //const surfaceScalarField& magSf = mesh.magSf();

        // Cell gradient of alpha
        if (nSmoothAlpha() > 0)
        {
            Info<<"Smoothing alpha"<<endl;
            volScalarField aveAlpha =  alpha1_;
            for (int nSmooth = 0; nSmooth < nSmoothAlpha(); nSmooth++)
            {
                aveAlpha = min(aveAlpha,0.95*aveAlpha + 0.05*fvc::average(aveAlpha));
            }
            gradAlphaPtr_.reset(new volVectorField(fvc::grad(aveAlpha, "nHat")));
        }
        else
        {
            gradAlphaPtr_.reset(new volVectorField(fvc::grad(alpha1_, "nHat")));
        }

        // Face gradient of alpha
        surfaceVectorField gradAlphaf(fvc::interpolate(gradAlphaPtr_()));

        // Face unit interface normal
        surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

        // Get time step size
        scalar deltaT(mesh.time().deltaT().value());

        // Get face-centered velocity
        surfaceVectorField Uf(fvc::interpolate(U_));

        // compute dynamic cAlpha based on Lee et. al 2015, 10th OFW
        // A Study of Computational Schemes for Six Degree-of-Freedom Motion of a Ship in Waves using OpenFOAM
        surfaceScalarField thetaf
        (
            IOobject
            (
                "thetaf",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        );

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        forAll(thetaf, faceI)
        {
            label own = owner[faceI];
            label nei = neighbour[faceI];

            vector d = C[nei] - C[own];

            thetaf[faceI] = ((d + (Uf[faceI]*deltaT)) & gradAlphaf[faceI])
                          / (mag(d + (Uf[faceI]*mag(deltaT))) * (mag(gradAlphaf[faceI]) + deltaN_.value()));
        }

        surfaceScalarField::Boundary& bthetaf = thetaf.boundaryFieldRef();
        surfaceVectorField::Boundary& bUf = Uf.boundaryFieldRef();
        surfaceVectorField::Boundary& bgradAlphaf = gradAlphaf.boundaryFieldRef();

        forAll(bthetaf, patchi)
        {
            scalarField& pthetaf = bthetaf[patchi];
            vectorField& pUf = bUf[patchi];
            vectorField& pgradAlphaf = bgradAlphaf[patchi];

            vectorField pd( bthetaf[patchi].patch().delta() );

            forAll(pthetaf, faceI)
            {
                pthetaf[faceI] = ((pd[faceI] + pUf[faceI]*deltaT) & pgradAlphaf[faceI])
                              / (mag(pd[faceI] + pUf[faceI]*deltaT) * (mag(pgradAlphaf[faceI]) + deltaN_.value()));
            }
        }

        //surfaceScalarField thetaf = ( (Sf + Uf*deltaT) & gradAlphaf )
        //                          / ( (magSf + mag(Uf)*deltaT) * (mag(gradAlphaf) + deltaN_) );

        cAlpha_ = Foam::sqrt(mag(cos(acos( thetaf ))));

        if (debug && mesh.time().write())
        {
            cAlpha_.write();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        IOobject
        (
            "cAlpha",
            alpha1.time().timeName(),
            alpha1.mesh()
        ),
        alpha1.mesh(),
        dimensionedScalar
        (
            "cAlpha",
            dimless,
            readScalar
            (
                alpha1.mesh().solution().solverDict(alpha1.name()).lookup("cAlpha")
            )
        )
    ),
    cAlphaModel_
    (
        alpha1.mesh().solution().solverDict(alpha1.name()).lookupOrDefault<word>
        ("cAlphaModel", "constant")
    ),
    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),
    alphaCoThreshold_
    (
        "alphaCoThreshold",
        alpha1.dimensions(),
        alpha1.mesh().time().controlDict().lookupOrDefault<scalar>("alphaCoThreshold", 0.05)
    ),
    nSmoothAlpha_
    (
        alpha1.mesh().solution().solverDict(alpha1.name()).lookupOrDefault<label>
        ("nSmoothAlpha", 0)
    ),
    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    read();
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterfaceCo() const
{
    return
        pos0(alpha1_ - alphaCoThreshold_)
       *pos0(1. - alphaCoThreshold_ - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
}


bool Foam::interfaceProperties::read()
{
    scalar cAlpha;
    alpha1_.mesh().solution().solverDict(alpha1_.name()).lookup("cAlpha") >> cAlpha;
    cAlpha_ = cAlpha;
    calculateCAlpha();
    sigmaPtr_->readDict(transportPropertiesDict_);

    alphaCoThreshold_ = dimensionedScalar
    (
        "alphaCoThreshold",
        alpha1_.dimensions(),
        alpha1_.mesh().time().controlDict().lookupOrDefault<scalar>("alphaCoThreshold", 0.05)
    );

    return true;
}


// ************************************************************************* //
