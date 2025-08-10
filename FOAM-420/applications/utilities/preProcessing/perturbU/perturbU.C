/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
// The FOAM Project // File: perturbU.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   Ucomponents
   \\  /           | Family: Utility
    \\/            |
    F ield         | FOAM version: 1.9.5
    O peration     |
    A and          | Copyright (C) 1991-1999 Henry G Weller
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
APPLICATION
    perturbU

DESCRIPTION
    initialise laminar velocity with superimposed streamwise streaks

AUTHOR
    Eugene de Villiers

-------------------------------------------------------------------------------
*/

#include "cfdTools/general/include/fvCFD.H"
#include "primitives/random/Random/Random.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "include/setRootCase.H"

#include "include/createTime.H"
#include "include/createMesh.H"

    Info<< "perturbU: generalised velocity perturbation implementation for "
         << "the initialisation of ducted well-resolved LES flows." << endl;

    //Xdir -> Ubar - streamwise
    //Ydir -> wallReflection vectors
    //Zdir -> cross product Xdir^Ydir

    //Ubar and Retau should be set in transportProperties
    //The choice of Retau is not critical as long as it is
    //approximately the right order of magnitude

    //A laminar background velocity profile is assumed
    //with maximum U at h = max(wall distance)
    //A laminar initial profile is essential since wall normal motion
    //of a mean turbulent profile without resolved turbulence will
    //diffuse the perturbations, preventing transition in some cases


    const wallDist& wd(wallDist::New(mesh));

    volScalarField yw(wd.y());
    const scalar h = max(yw.primitiveField());

    //local yDir
    const volVectorField yDir = wd.n();


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U" << endl;
    volVectorField U(Uheader, mesh);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        "nu", transportProperties.lookup("nu")
    );
    dimensionedVector Ubar
    (
        "Ubar", transportProperties.lookup("Ubar")
    );
    dimensionedScalar Retau
    (
        "Retau", transportProperties.lookup("Retau")
    );

    vector xDir = Ubar.value()/mag(Ubar.value());

    Switch pipe = transportProperties.lookupOrDefault<Switch>("pipe", true);


    Info<< "Re(tau) = " << Retau << endl;
    const scalar utau = Retau.value()*nu.value()/h;
    Info<< " u(tau) = " << utau << endl;


    //wall normal circulation
    const scalar duplus = Ubar.value().x()*0.25/utau;
    //spanwise wavenumber: spacing z+ = 200
    const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/200.0);
    const scalar sigma = 0.00055;
    //streamwise wave number: spacing x+ = 500
    const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/500.0);
    const scalar epsilon = Ubar.value().x()/200.0;

    const vectorField& centres(mesh.C());

    //pipe component
    vector Cgeom
        = gSum(mesh.C().field()*mesh.V().field())/gSum(mesh.V().field());
    vector e2g(0, 1, 0);

    if (pipe)
    {
        if ((xDir & e2g) > (1 - SMALL))
        {
            e2g = vector(0, 0, 1);
        }

        e2g -= (e2g & xDir) * xDir;

        e2g /= mag(e2g);

    }

    Random perturbation(1234567);

    forAll(centres, celli)
    {
        //add a small random component to enhance symmetry breaking
        scalar deviation=1.0 + 0.1*perturbation.GaussNormal<scalar>();
        const vector& cCentre = centres[celli];

        vector zDir = xDir^yDir[celli];
        zDir /= mag(zDir);

        scalar zplus = (cCentre & zDir)*Retau.value()/h;
        scalar yplus = yw[celli]*Retau.value()/h;
        scalar xplus = (cCentre & xDir)*Retau.value()/h;

        // laminar parabolic profile
        if (pipe) //circular duct
        {
            U[celli] = 2.0*Ubar.value()*(1 - sqr(1 - yw[celli]/h));

            vector r = cCentre - Cgeom;
            r -= (r & xDir) * xDir;
            scalar magr = max(SMALL, mag(r));

            scalar theta = ::acos((r & e2g)/magr);

            zplus = theta * magr * Retau.value()/h;
        }
        else // parallel flat plates
        {
            U[celli] = 3.0*Ubar.value() * (yw[celli]/h - 0.5*sqr(yw[celli]/h));
        }

        // streak streamwise velocity
        U[celli] +=
            xDir*(utau * duplus/2.0) * (yplus/40.0)
            * Foam::exp(-sigma * Foam::sqr(yplus) + 0.5)
            * Foam::cos(betaPlus*zplus)*deviation;

        // streak spanwise perturbation
        U[celli] +=
            zDir * epsilon
          * Foam::sin(alphaPlus*xplus)
          * yplus
          * Foam::exp(-sigma*Foam::sqr(yplus))
          * deviation;

    }

    #include "cfdTools/incompressible/createPhi.H"

    U.correctBoundaryConditions();

    U.write();

    Info<< endl;

    return(0);
}


// ************************************************************************* //
