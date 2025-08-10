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
    (c) held by original author
    (c) 2022 Esi Ltd.

Class
    interface

SourceFiles
    interface.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de> (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "interface.H"
#include "alphaContactAngle/alphaContactAngle/alphaContactAngleFvPatchScalarField.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvc/fvcSnGrad.H"


#include<sys/types.h>
#include<iostream>
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interface::convertToRad =
                Foam::constant::mathematical::pi / 180.0;

// * * * * * * * * * * * * * * * Private Member Function  * * * * * * * * * * //

//- Creates a volScalarField isInterface which is initialized with
//  calculateIsInterface(volScalarField&) to be called in the interface constructor
void Foam::interface::initIsInterface()
{
    calculateIsInterface(isInterface_);
}

//- Marks cells which will be handled as interface cells. Two Models are
//  selectable one with the raw alpha value and another using its gradient
//  For better control of the marked cell band-width additional cell layers
//  can be added
void Foam::interface::calculateIsInterface(Foam::volScalarField& isInterface)
{
    if (isInterfaceMethod_ == "snGradAlpha")
    {
        isInterface = 0.0;
        isInterface.primitiveFieldRef()
            =   pos0(
                    fvc::average(
                        Foam::mag(Foam::fvc::snGrad(alpha1_))
                        / alpha1_.mesh().deltaCoeffs()
                    )
                    - isInterfaceThreshold_
                );

    }
    else if (isInterfaceMethod_ == "alpha")
    {
        isInterface = 0.0;
        forAll(alpha1_, iCell)
        {
          if (alpha1_[iCell] > isInterfaceThreshold_ && alpha1_[iCell] < (1- isInterfaceThreshold_))
          {
            isInterface[iCell] = 1;
          }
        }
    }
    else
    {
        Info<< "In interface::calculateIsInterface()\n"
             << "choose an appropriate isInterfaceMethod\n"
             << "e.g. snGradAlpha" << endl;
    }

    //- Get a wider and smoother stencil
    for (int i = 0; i < isInterfaceAddN_; i++)
    {
        isInterface =
            pos0(fvc::average(fvc::interpolate(isInterface)) - SMALL);
    }

    //- set internal isInterfaceMarkers on boundary field
    forAll(isInterface.boundaryField(), iPatch)
    {
        scalarField& isIp = isInterface.boundaryFieldRef()[iPatch];
        isIp = isInterface.boundaryField()[iPatch].patchInternalField();
    }
}


void Foam::interface::calculateCAlpha()
{
}


void Foam::interface::initialise()
{
    read();

    if (transpProp_.found("reconstruct"))
    {
        reconI_ =
        (
           reconstruct::New
           (
               transpProp_.subDict("reconstruct").lookup("reconstructModel"),
               alpha1_,
               transpProp_,
               isWallPatch_,
               isInterface_
           )
        );
    }
    else
    {
        WarningInFunction
            << "transportProperties subDictionary \"reconstruct\" not found. "
            << "Please double-check your setup."
            << nl << endl;
    }

    if (transpProp_.found("curvature"))
    {
        K_ =
        (
           curvature::New
           (
               transpProp_.subDict("curvature").lookup("curvatureModel"),
               alpha1_,
               reconI_->nHatv(),
               reconI_->nHatfv(),
               reconI_->interfaceDensity(),
               isWallPatch_,
               transpProp_,
               isInterface_
           )
        );
    }
    else
    {
        WarningInFunction
            << "transportProperties subDictionary \"curvature\" not found. "
            << "Please double-check your setup."
            << nl << endl;
    }

    initIsInterface();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interface::interface
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
    :
    transpProp_(dict),

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

    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0), transpProp_.lookup("sigma")),

    deltaN_
    (
       "deltaN",
       1e-8 / pow(average(alpha1.mesh().V()), 1.0 / 3.0)
    ),
    alphaCoThreshold_
    (
        "alphaCoThreshold",
        alpha1.dimensions(),
        alpha1.mesh().time().controlDict().lookupOrDefault<scalar>("alphaCoThreshold", 0.05)
    ),

    alpha1_(alpha1),
    U_(U),

    isWallPatch_(setIsWallPatch(alpha1)),

    isInterfaceMethod_(),

    isInterfaceThreshold_(),

    isInterfaceAddN_(),

    isInterface_
    (
        IOobject
        (
           "isInterface",
           alpha1_.mesh().time().timeName(),
           alpha1_.mesh(),
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       alpha1_.mesh(),
       dimensionedScalar("isInterface", dimless, 0)
    ),

    reconI_(),

    K_()
{
    initialise();
    Info<< "interface::interface() Interface is created " << endl;
}


//- Determines which patches are alphacContactAngle patches
Foam::List<bool> Foam::interface::setIsWallPatch(const volScalarField& alpha) const
{
    Foam::List<bool> isWall(1, false);

    //- count number of patches
    scalar count = 0;
    forAll(alpha.boundaryField(), iPatch)
    {
        count++;
    }
    isWall.resize(count, false);

    //- true if patch is a wall
    forAll(alpha1_.boundaryField(), iPatch)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[iPatch]))
        {
            Info<< "interface::setIsWallPatch  yes " << iPatch <<"  "<< alpha1_.boundaryField()[iPatch].type() << endl;
            isWall[iPatch] = true;
        } else {
            Info<< "interface::setIsWallPatch  no  " << iPatch <<"  "<< alpha1_.boundaryField()[iPatch].type() << endl;
        }
    }

    return isWall;
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interface::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
)
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

            scalarField a((b1 - a12*b2)/(det + VSMALL));
            scalarField b((b2 - a12*b1)/(det + VSMALL));

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
        else if (
            (!isWallPatch_[patchi])
            && (alpha1_.boundaryField()[patchi].type() != "processor")
            && (alpha1_.boundaryField()[patchi].type() != "cyclic")
        )
        {
            //Info<<"new correction for patch: " << boundary[patchi].name() << "\n"
            //   << " and type " << boundary[patchi].type() << endl;
            fvsPatchVectorField& nHatp = nHatb[patchi];

            //- set contact angle to 90°
            vectorField nf( alpha1_.mesh().boundary()[patchi].nf() );
            nHatp -= (nHatp & nf)*nf;
            nHatp /= (mag(nHatp) + VSMALL);

            //- NOT Necessary if complexDivCorrection is selected for curvature calculation:
            reconI_->nHatv().boundaryFieldRef()[patchi] = nHatp;
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::interface::surfaceTensionForce() const
{
    return K_->surfaceTensionForce(sigma_);
}

Foam::tmp<Foam::volScalarField>
Foam::interface::nearInterface() const
{
    return isInterface_;
}

Foam::tmp<Foam::volScalarField>
Foam::interface::nearInterfaceCo() const
{
    return
        pos0(alpha1_ - alphaCoThreshold_)
       *pos0(1. - alphaCoThreshold_ - alpha1_);
}

bool Foam::interface::read()
{
    scalar cAlpha;
    alpha1_.mesh().solution().solverDict(alpha1_.name()).lookup("cAlpha") >> cAlpha;
    cAlpha_ = cAlpha;
    calculateCAlpha();

    alphaCoThreshold_ = dimensionedScalar
    (
        "alphaCoThreshold",
        alpha1_.dimensions(),
        alpha1_.mesh().time().controlDict().lookupOrDefault<scalar>("alphaCoThreshold", 0.05)
    );

    if (transpProp_.found("interface"))
    {
        const dictionary& dict = transpProp_.subDict("interface");
        isInterfaceMethod_ = dict.lookupOrDefault<word>("isInterfaceMethod", "snGradAlpha");
        isInterfaceThreshold_= dict.lookupOrDefault<scalar>("isInterfaceThreshold", 0.05);
        isInterfaceAddN_ = dict.lookupOrDefault<scalar>("isInterfaceAddN", 1);
    }
    else
    {
        WarningInFunction
            << "transportProperties subDictionary \"interface\" not found. "
            << "Please double-check your setup."
            << nl << endl;
    }

    return true;
}

// ************************************************************************* //
