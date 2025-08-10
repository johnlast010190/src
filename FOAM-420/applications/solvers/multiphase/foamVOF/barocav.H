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
    (c) 2019-2020 Esi Ltd.

Application
    foamVOF

Group
    grpMultiphaseSolvers

Description
    Solver for multi-phase flows with two or more fluids. This file contains
    implementation details of the barotropic cavitation model

\*---------------------------------------------------------------------------*/

#include "vof.h"
#include "barocav.h"

using namespace Foam;

BaroCav::BaroCav
(
    const Time& time,
    const dynamicFvMesh& mesh,
    const MatProp &prop
)
:
FoamVOF
(
    time,
    mesh,
    prop
),
maxAcousticCo(1.0)
{
    name_="baroCavitate";
    correctPhi_=false;
}

volScalarField& BaroCav::alpha1()
{
    return mixture().alpha1();
}

volScalarField& BaroCav::alpha2()
{
    return mixture().alpha2();
}

void BaroCav::createFields(pimpleControl &pimple)
{
    Info<< "Reading thermodynamicProperties\n" << endl;

//    const Time &runTime=runTime();

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime().constant(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar psil0
    (
        "psil",
        dimCompressibility,
        thermodynamicProperties
    );
    psil.dimensions().reset(dimensionSet(dimCompressibility));
    psil=psil0;

    dimensionedScalar rholSat0
    (
        "rholSat",
        dimDensity,
        thermodynamicProperties
    );
    rholSat.dimensions().reset(dimensionSet(dimDensity));
    rholSat=rholSat0;

    dimensionedScalar psiv0
    (
        "psiv",
        dimCompressibility,
        thermodynamicProperties
    );
    psiv.dimensions().reset(dimensionSet(dimCompressibility));
    psiv=psiv0;

    dimensionedScalar pSat0
    (
        "pSat",
        dimPressure,
        thermodynamicProperties
    );
    pSat.dimensions().reset((dimPressure));
    pSat=pSat0;

    dimensionedScalar rhovSat0("rhovSat", psiv*pSat);

    rhovSat.dimensions().reset(rhovSat0.dimensions());
    rhovSat=rhovSat0;

    dimensionedScalar rhol00("rhol0", rholSat - pSat*psil);
    rhol0.dimensions().reset(rhol00.dimensions());
    rhol0=rhol00;

    dimensionedScalar rhoMin0
    (
        "rhoMin",
        dimDensity,
        thermodynamicProperties
    );

    rhoMin.dimensions().reset(dimDensity);
    rhoMin=rhoMin0;

    Info<< "Reading field p\n" << endl;
    autoPtr<volScalarField> pPtr
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
    pPtr.ptr()->store();

    autoPtr<volScalarField> rhoPtr
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
    rhoPtr.ptr()->store();

    Info<< "Reading field U\n" << endl;
    autoPtr<volVectorField> UPtr
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    UPtr.ptr()->store();

if (!isA<staticFvMesh>(mesh()))
{
    Info<< "Reading/calculating  field Uf\n" << endl;
    autoPtr<surfaceVectorField> UfPtr
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U())
        )
    );
    UfPtr.ptr()->store();

// createPcorrTypes.H

   pcorrTypes=
   new wordList
    (
        p().boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p().boundaryField().size(); i++)
    {
        if (p().boundaryField()[i].fixesValue())
        {
            pcorrTypes()[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

} //dynami mesh

    //#include "cfdTools/incompressible/createPhi.H"
    Info<< "Reading/calculating face flux field phi\n" << endl;

    autoPtr<surfaceScalarField> phiPtr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(U())
        )
    );

    phiPtr.ptr()->store();
    mesh().schemes().setFluxRequired(p().name());

    // Mass flux (corrected by rhoEqn.H)
    autoPtr<surfaceScalarField> rhoPhiPtr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rho())*phi()
        )
    );
    rhoPhiPtr.ptr()->store();

    Info<< "Reading transportProperties\n" << endl;

    mixture2.reset(new incompressibleTwoPhaseMixture(U(), phi()));



    //volScalarField& alphav(mixture().alpha1());
    alphav().oldTime();

    //volScalarField& alphal(mixture().alpha2());

    Info<< "Creating compressibilityModel\n" << endl;

    //autoPtr<barotropicCompressibilityModel> psiModel =
    psiModel =
    barotropicCompressibilityModel::New
    (
        thermodynamicProperties,
        alphav()
    );


//    const volScalarField& psi = psiModel->psi();


    rho().forceAssign(max
    (
        psi()*p()
    + alphal()*rhol0
    + ((alphav()*psiv + alphal()*psil) - psi())*pSat,
        rhoMin
    ));

    mesh().schemes().setFluxRequired(p().name());
    mesh().schemes().setFluxRequired(rho().name());

    createTfields(pimple);
    return;
}

void BaroCav::CourantNo()
{
    CoNum = 0.0;
    meanCoNum = 0.0;
    acousticCoNum = 0.0;


    if (mesh().nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi()))().primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/mesh().V().field())*dt();

        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh().V().field()))*dt();

        acousticCoNum = 0.5*gMax
        (
            fvc::surfaceSum
            (
                fvc::interpolate(scalar(1)/sqrt(psi()))*mesh().magSf()
            )().primitiveField()/mesh().V().field()
        )*dt();
    }

    Info<< "phi Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " acoustic max: " << acousticCoNum
        << endl;

    return;
}

void BaroCav::UEqnSolve
(
    fvVectorMatrix& UEqn,
    pimpleControl &pimple,
    incompressible::turbulenceModel& turbulence
)
{
    UEqn =fvVectorMatrix
    (
        fvm::ddt(rho(), U())
      + fvm::div(rhoPhi(), U())
      + turbulence.divDevRhoReff(rho(), U())
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p()));
    }

    Info<< "max(U) " << max(mag(U())).value() << endl;
    return;
}


void BaroCav::rhoEqnSolve()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho())
      + fvm::div(phi(), rho())
    );

    rhoEqn.solve();

    rhoPhi() = rhoEqn.flux();

    Info<< "max-min rho: " << max(rho()).value()
        << " " << min(rho()).value() << endl;

    rho().forceAssign(max(rho(), rhoMin));
    return;
}

void BaroCav::correctAlphav()
{
    alphav() =
        max
        (
            min
            (
                (rho() - rholSat)/(rhovSat - rholSat),
                scalar(1)
            ),
            scalar(0)
        );

    alphal() = 1.0 - alphav();

    Info<< "max-min alphav: " << max(alphav()).value()
        << " " << min(alphav()).value() << endl;

    psiModel->correct();
    return;
}

void BaroCav::pEqnSolve
(
    fvVectorMatrix& UEqn,
    pimpleControl &pimple,
    incompressible::turbulenceModel& turbulence
)
{
    if (pimple.nCorrPIMPLE() == 1)
    {
        p() =
        (
            rho()
          - alphal()*rhol0
          - ((alphav()*psiv + alphal()*psil) - psi())*pSat
        )/psi();
    }

    surfaceScalarField rhof("rhof", fvc::interpolate(rho()));

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho()*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p()));
    if (isA<staticFvMesh>(mesh()))
    {
        phi() = fvc::flux(HbyA)
            + rhorAUf*fvc::ddtCorr(U(), phi());
    }
    else
    {
        phi() = fvc::flux(HbyA)
         + rhorAUf*fvc::ddtCorr(U(), Uf());
        fvc::makeRelative(phi(), U());
    }

    surfaceScalarField phiGradp(rhorAUf*mesh().magSf()*fvc::snGrad(p()));

    phi() -= phiGradp/rhof;

    while (pimple.correctNonOrthogonal())
    {
        if (isA<staticFvMesh>(mesh()))
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi(), p())
                - (rhol0 + (psil - psiv)*pSat)*fvc::ddt(alphav()) - pSat*fvc::ddt(psi())
                + fvc::div(phi(), rho())
                + fvc::div(phiGradp)
                - fvm::laplacian(rhorAUf, p())
            );

            pEqn.solve(mesh().solution().
                solver(p().select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi() += (phiGradp + pEqn.flux())/rhof;
            }
        }
        else
        {
            fvScalarMatrix pEqn
            (
                fvc::ddt(rho())
                + psi()*correction(fvm::ddt(p()))
                + fvc::div(phi(), rho())
                + fvc::div(phiGradp)
                - fvm::laplacian(rhorAUf, p())
            );

            pEqn.solve(mesh().solution().
                solver(p().select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi() += (phiGradp + pEqn.flux())/rhof;
            }
        }
    }

    Info<< "Predicted p max-min : " << max(p()).value()
        << " " << min(p()).value() << endl;

    if (isA<staticFvMesh>(mesh()))
    {
        rho().forceAssign(max
        (
            psi()*p()
            + alphal()*rhol0
            + ((alphav()*psiv + alphal()*psil) - psi())*pSat,
            rhoMin
        ));
    }
    else
    {
        volScalarField rho0(rho() - psi()*p());
        rho().forceAssign(max(rho0 + psi()*p(), rhoMin));
    }

    correctAlphav();

    p() =
    (
        rho()
      - alphal()*rhol0
      - ((alphav()*psiv + alphal()*psil) - psi())*pSat
    )/psi();

    p().correctBoundaryConditions();

    Info<< "Phase-change corrected p max-min : " << max(p()).value()
        << " " << min(p()).value() << endl;

    // Correct velocity

    U() = HbyA - rAU*fvc::grad(p());

    // Remove the swirl component of velocity for "wedge" cases
    if (pimple.dict().found("removeSwirl"))
    {
        label swirlCmpt(readLabel(pimple.dict().lookup("removeSwirl")));

        Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;
        U().field().replace(swirlCmpt, 0.0);
    }

    U().correctBoundaryConditions();

    Info<< "max(U) " << max(mag(U())).value() << endl;

    if (!isA<staticFvMesh>(mesh()))
    {
        Uf() = fvc::interpolate(U());
        surfaceVectorField n(mesh().Sf()/mesh().magSf());
        Uf() += n*(phi()/mesh().magSf() - (n & Uf()));
    }

    return;
}

void BaroCav::correctPhi(pimpleControl &pimple)
{
    //FoamVOF::correctPhi(pimple);
    correctUphiBCs(U(), phi());

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pcorr", p().dimensions(), 0.0),
        pcorrTypes()
    );

    surfaceScalarField rhof(fvc::interpolate(rho(), "div(phi,rho)"));
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    mesh().schemes().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::ddt(rho()) + fvc::div(phi()*rhof)
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi() -= pcorrEqn.flux()/rhof;
        }
    }

    return;
}


void BaroCav::updateFields
(
    pimpleControl &pimple,
    incompressible::turbulenceModel &turbulence,
    fv::options& fvOptions
)
{
    readTimeControls();
    maxAcousticCo= readScalar(runTime().controlDict().lookup("maxAcousticCo"));
    correctPhi_ = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

    CourantNo();

    if (adjustTimeStep)
    {
      scalar maxDeltaTFact =
        min(maxCo/(CoNum + SMALL), maxAcousticCo/(acousticCoNum + SMALL));

      scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

      const_cast<Time&> (runTime()).setDeltaT
      (
        min
        (
            deltaTFact*runTime().deltaTValue(),
            maxDeltaT
        )
      );

      Info<< "deltaT = " <<  runTime().deltaTValue() << endl;

      //exit Foam if deltaT is below user secified minDeltaT
      if (runTime().deltaTValue() < minDeltaT)
      {
          FatalErrorInFunction
              << "minDeltaT = " << minDeltaT << nl
              << "The computed deltaT is lower than the minDeltaT"
              << exit(FatalError);
      }

   }


    const_cast<Time&> (runTime_)++;
    Info<< "Time = " << runTime().timeName() << nl << endl;
    if (!isA<staticFvMesh>(mesh()))
    {
        // Do any mesh changes
        const_cast<dynamicFvMesh&> (mesh()).update();

        if (mesh().changing() && correctPhi())
        {
            // Calculate absolute flux from the mapped surface velocity
            phi() = mesh().Sf() & Uf();

            //#include "correctPhi.H"
            correctPhi(pimple);
            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi(), U());
        }
    }
    // --- Pressure-velocity PIMPLE corrector loop

    while (pimple.loop())
    {
        rhoEqnSolve();
        correctAlphav();

        fvVectorMatrix UEqn(U(), dimensionSet(1, 1, -2, 0, 0, 0, 0));
        UEqnSolve(UEqn,pimple,turbulence);

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            pEqnSolve(UEqn,pimple,turbulence);
        }

        if (pimple.turbCorr())
        {
            turbulence.correct();
        }
    }
    return;
}
