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

Class
    reconstruct

SourceFiles
    reconstruct.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Vignesh Thammanna Gurumurthy < vignesh@sla.tu-darmstadt.de >

    All rights reserved.


S-SCLSVOF - Density Balanced - Method Developement:

    Albadawi A, Donoghue DB, Robinson AJ, Murray DB, Delauré YMC.
    Influence of surface tension implementation in volume of fluid and coupled volume of fluid with level set methods for bubble growth and detachment. International
    Journal of Multiphase Flow 2013; 53:11–28.

    T. YAMAMOTO, Y. OKANO AND S. DOST
    Validation of the S-CLSVOF method with the density-scaled balanced continuum surface force model in multiphase systems coupled with thermocapillary flows
    Int. J. Numer. Meth. Fluids 2017; 83:223-244.

   Yokoi K.
    A density-scaled continuum surface force model within a balanced force formulation. Journal of
    Computational Physics 2014; 278:221–228.

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de> (main developer).

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


#include "reconstruct/reconstructLevelSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug( reconstructLevelSet, 0 );
addToRunTimeSelectionTable( reconstruct, reconstructLevelSet, dictionary );


void reconstructLevelSet::calculateDelta()
{
    forAll(mesh_.cells(),celli)
    {
        if (mag(psi_[celli]) > epsilon_[celli])
        {
            delta_[celli] = 0.0;
        }
        else
        {
            delta_[celli] = 1.0/(2.0*epsilon_[celli])*(1.0+cos(pi*psi_[celli]/epsilon_[celli]));
        }
    }
    delta_.correctBoundaryConditions();

    if (scaled_)
    {
        delta_ = 2.0*H_*delta_;
        delta_.correctBoundaryConditions();
    }
}

void reconstructLevelSet::calculateH()
{
    forAll(mesh_.cells(),celli)
    {
        if (psi_[celli] < -epsilon_[celli])
        {
            H_[celli] = 0.0;
        }
        else if (epsilon_[celli] < psi_[celli])
        {
            H_[celli] = 1.0;
        }
        else
        {
            H_[celli] =   0.5 *(  1
                                  + psi_[celli]/epsilon_[celli]
                                  + sin(pi*psi_[celli]/epsilon_[celli])/pi
                                );
        }
    }
    H_.correctBoundaryConditions();

    if (scaled_)
    {
        forAll(mesh_.cells(),celli)
        {
            if (psi_[celli] < -epsilon_[celli])
            {
                Hscaled_[celli] = 0.0;
            }
            else if (epsilon_[celli] < psi_[celli])
            {
                Hscaled_[celli] = 1.0;
            }
            else
            {
                Hscaled_[celli] = 0.5
                                  * (  0.5
                                     + psi_[celli]/epsilon_[celli]
                                     + psi_[celli]*psi_[celli]/(2.0*epsilon_[celli]*epsilon_[celli])
                                     - (cos(2.0*pi*psi_[celli]/epsilon_[celli])-1)/(4.0*pi*pi)
                                     +   sin(pi*psi_[celli]/epsilon_[celli])
                                       * (epsilon_[celli]+psi_[celli])/(pi*epsilon_[celli])
                                    );
            }
        }
        Hscaled_.correctBoundaryConditions();
    }
}

void reconstructLevelSet::initConstants()
{
    deltaX_.primitiveFieldRef() = pow(mesh_.V(), 1.0/3.0);

    //- use minimal cell lenght as global deltaX to avoid high occurences
    //  of delta and H because of different cell refinement levels.
    //  Should not influence the interface curvature if the whole
    //  interface shares the same refinement level
    scalar minDelta = gMin(deltaX_);
    gamma_ = 0.75*minDelta;
    deltaTau_ = 0.1*minDelta;
    epsilon_= epsC_*minDelta;

}

void reconstructLevelSet::solveReinitializationEqn()
{
    // solve Level-Set function as the re-initialization equation
    Info<< "solve the reinitialization equation" << nl << endl;

    //- calculate psi0_
    psi0_ == (2.0*alpha1_ - 1.0)*gamma_;
    psi0_.correctBoundaryConditions();

    //- correct psi_
    psi_ == psi0_;

    dimensionedScalar ds("ds", dimensionSet(0,1,0,0,0,0,0), 1.0);
    for (int corr=0; corr<epsC_*10; corr++)
    {
      psi_ == psi_ + psi0_/mag(psi0_+VSMALL)*(1.0-mag(fvc::grad(psi_)*ds))*deltaTau_;
      psi_.correctBoundaryConditions();
    }

    //- calculate H and if selected Hscaled
    calculateH();

    //- calculate delta which if selected represents deltaScaled
    calculateDelta();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reconstructLevelSet::reconstructLevelSet
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& transpProp,
    const List<bool>& isWallPatch,
    const volScalarField& isInterface
)
:
    reconstruct( name, alpha, transpProp, isWallPatch, isInterface ),

    mesh_(alpha.mesh()),

    psi_
    (
        IOobject
        (
            "psi",
            alpha.mesh().time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh()
    ),

    gradPsif_
    (
        IOobject
        (
            "gradPsif",
            alpha.mesh().time().timeName(),
            alpha.mesh()
        ),
        alpha.mesh(),
        dimensionedVector("gradPsif", dimensionSet(0,-1,0,0,0,0,0), vector(0,0,0))
    ),

    psi0_
    (
        IOobject
        (
            "psi0",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),
        dimensionedScalar("psi0", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
        //psi_.boundaryField().types()
    ),

    delta_
    (
        IOobject
        (
            "delta",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        psi_.mesh(),
        dimensionedScalar("delta", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    H_
    (
        IOobject
        (
            "H",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),
        dimensionedScalar("H", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
        //psi_.boundaryField().types()
    ),

    Hscaled_
    (
        IOobject
        (
            "Hscaled",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),
        dimensionedScalar("Hscaled", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
        //psi_.boundaryField().types()
    ),

    deltaX_
    (
        IOobject
        (
            "deltaX",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),//, pow(alpha.mesh().V(), 1.0/3.0)
        dimensionedScalar("deltaX", dimensionSet(0,0,0,0,0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    deltaTau_
    (
        IOobject
        (
            "deltaTau_",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),//,
        dimensionedScalar("deltaTau_", dimensionSet(0,0,0,0,0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    gamma_
    (
        IOobject
        (
            "gamma",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(), //,        0.75*deltaX_
        dimensionedScalar("Hscaled", dimensionSet(0,0,0,0,0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    epsC_(readScalar(transpProp.subDict("reconstruct").lookup("epsC"))),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            alpha.mesh().time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(), //,         epsC_*deltaTau_
        dimensionedScalar("Hscaled", dimensionSet(0,0,0,0,0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    scaled_(transpProp.subDict("reconstruct").lookup("scaled"))

{
    //- initialize LevelSet constants
    initConstants();

    //- initialize nHatfv
    reconstructInterface();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reconstructLevelSet::reconstructInterface()
{
    if (mesh_.changing())
    {
        initConstants();
    }

    solveReinitializationEqn();

    // Cell gradient of psi
    const volVectorField gradPsi(fvc::grad(psi_));

    // Interpolated face-gradient of psi
    gradPsif_ = fvc::interpolate(gradPsi);

    //- calculate the face unit interface normal field
    nHatfv_ = gradPsif_/(mag(gradPsif_) + deltaN_);
    nHatv_  = gradPsi/(mag(gradPsi) + deltaN_);

    //- calculate interface density
    interfaceDensity_ = mag(gradPsi); //TODO is this correct?
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
