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
    (c) 2019-2022 Esi Ltd.

Application
    foamVOF

Group
    grpMultiphaseSolvers

Description
    Solver for multi-phase flows with two or more fluids. This file contains
    implementation details of the member functions of the FoamVOF class

\*---------------------------------------------------------------------------*/

#include "vof.h"
#include "cavitate.h"

using namespace Foam;

FoamVOF::FoamVOF
(
    const Time& time,
    const dynamicFvMesh& mesh,
    const MatProp &prop0
)
:
    mesh_(mesh),
    runTime_(time),
    prop_(prop0),
    mfluid_(runTime_.controlDict().lookupOrDefault<label>("mfluid", 2)),
    solveT_(runTime_.controlDict().lookupOrDefault<label>("solveT", 0)),
    baroCavitate_(runTime_.controlDict().lookupOrDefault<label>("baroCavitate", 0)),
    algorithm_(runTime_.controlDict().lookupOrDefault<label>("algorithm", PIMPLE)),
   // use_rgh_(runTime_.controlDict().lookupOrDefault<label>("use_rho_gh", true)),
    cavmodel_(runTime_.controlDict().lookupOrDefault<label>("cavitation", 0)),
    moveMeshOuterCorrectors(false)
    {
        name_="foamVOF";
        mules_=true;
        g0_=0.0;
        if (!useRhogh())
        {
            g0_=9.81;
        }
        use_rgh_=runTime_.controlDict().lookupOrDefault<label>("use_rho_gh", 1);
    }

void FoamVOF::createModels()
{
    VOFModels *mdl= nullptr;
    if (cavitation()>0)
    {
      Info<<"Create cavitation model..."<<endl;
      mdl=new Cavitate(prop_,*this);
      mules_=false;
    }
    if (mdl)
    {
        models.append(mdl);
    }
    Info<<"Number of physical models created:"<<models.size()<<endl;

    return;
}

bool FoamVOF::hasModel()
{
    if (models.size()>0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


surfaceScalarField& FoamVOF::alphaPhi()
{
    return  alphaPhi_.ref();
}


volScalarField& FoamVOF::alpha1()
{
    return  mixture->alpha1();
}

volScalarField& FoamVOF::alpha2()
{
    return  mixture->alpha2();
}

const dimensionedScalar& FoamVOF::rho1()
{
    return mixture->rho1();
}

const dimensionedScalar& FoamVOF::rho2()
{
    return mixture->rho2();
}

void FoamVOF::correctAlphav()
{
    return;
}

void FoamVOF::createAlphaFluxes()
{
    IOobject alphaPhiHeader
  (
    "alphaPhi",
    runTime().timeName(),
    mesh(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  );

 alphaRestart = alphaPhiHeader.typeHeaderOk<surfaceScalarField>(true);

 // MULES flux from previous time-step
 alphaPhi_=tmp<surfaceScalarField>
 (
    new surfaceScalarField
    (
     alphaPhiHeader,
     phi()*fvc::interpolate(alpha1())
    )
 );

    return;
}

void FoamVOF::createTfields
(
    pimpleControl &pimple
)
{
    if (solveT()==0)
    {
      return;
    }

    Info<< "Reading field T\n" << endl;
    autoPtr<volScalarField> TPtr
    (
        new volScalarField
        (
        IOobject
        (
        "T",
        runTime().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
        mesh()
        )
    );

    TPtr.ptr()->store();

    // thermal capacity cp
    dimensionedScalar cp1("cp", dimEnergy/dimTemperature/dimMass ,  prop().cp1());
    dimensionedScalar cp2("cp", dimEnergy/dimTemperature/dimMass ,  prop().cp2());
    autoPtr<volScalarField> cpPtr
    (
    new volScalarField
    (
        IOobject
        (
        "cp",
        runTime().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        alpha1()*cp1+alpha2()*cp2
        )
        );

    cpPtr.ptr()->store();
    Info<<"finished reading Cp"<<endl;
    //thermal conductivity

    dimensionedScalar K1("K",  dimensionSet(1,1,-3,-1,0,0,0) , prop().k1());
    dimensionedScalar K2("K",  dimensionSet(1,1,-3,-1,0,0,0) , prop().k2());

    autoPtr<volScalarField> KPtr
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),

            alpha1()*K1+alpha2()*K2
        )

  );

   KPtr.ptr()->store();
   return;
 }

void FoamVOF::createFields(pimpleControl &pimple)
{

    LTS = fv::localEulerDdt::enabled(mesh());
    if (LTS)
    {
        Info<< "Using LTS" << endl;
        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("one", dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }

    Info<< "Reading field p_rgh\n" << endl;

    autoPtr<volScalarField> p_rghPtr
    (
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    p_rghPtr.ptr()->store();

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
    }

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

    Info<< "Reading transportProperties\n" << endl;

    mixture.reset(new immiscibleIncompressibleTwoPhaseMixture(U(), phi()));

    volScalarField& alpha1(mixture->alpha1());
    volScalarField& alpha2(mixture->alpha2());

    const dimensionedScalar& rho1 = mixture->rho1();
    const dimensionedScalar& rho2 = mixture->rho2();

    // Need to store rho for ddt(rho, U)
    autoPtr<volScalarField> rhoPtr
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT
            ),
            alpha1*rho1 + alpha2*rho2
        )
    );

    rhoPtr.ptr()->store();
    rho().oldTime();
    createTfields(pimple);

    //volume transfer rate due to interphase mass transfer
    if (cavitation())
    {
        dimensionedScalar su0("su",  dimensionSet(0,0,-1,0,0,0,0) , 0.0);
        dimensionedScalar sp0("sp",  dimensionSet(0,0,-1,0,0,0,0)/dimPressure , 0.0);

        autoPtr<volScalarField> volTrRatePtr
        (
            new volScalarField
            (
                IOobject
                (
                    "volTrRate",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                alpha1*su0
            )
        );

        volTrRatePtr.ptr()->store();

        autoPtr<volScalarField> vdotPPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "vdotP",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                alpha1*sp0
            )
        );

        vdotPPtr.ptr()->store();
    }

    // Mass flux
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

    Info<< "\nReading g" << endl;

    autoPtr<uniformDimensionedVectorField> gPtr
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime().constant(),
                mesh(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    gPtr.ptr()->store();

    Info<< "\nReading hRef" << endl;
    autoPtr<uniformDimensionedScalarField> hRefPtr
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime().constant(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, 0)
        )
    );

    hRefPtr.ptr()->store();

    Info<< "Calculating field g.h\n" << endl;

    dimensionedScalar ghRef
    (
        mag(g().value()) > SMALL
        ? g() & (cmptMag(g().value())/mag(g().value()))*hRef()
        : dimensionedScalar("ghRef", g().dimensions()*dimLength, 0)
    );

    ghRefv=ghRef.value();


    autoPtr<volScalarField> ghPtr
    (
        new volScalarField
        (
            IOobject
            (
                "gh",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (g() & mesh().C()) - ghRef
        )
    );
    ghPtr.ptr()->store();


    autoPtr<surfaceScalarField> ghfPtr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ghf",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (g() & mesh().Cf()) - ghRef
        )
    );
    ghfPtr.ptr()->store();

    if (mixture->found("frameAcceleration"))
    {
        frameAcceleration.reset
        (
            Function1<vector>::New("frameAcceleration", mixture()).ptr()
        );
    }


    autoPtr<volScalarField> pPtr
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p_rgh() + rho()*gh()
        )
    );
    pPtr.ptr()->store();

    pRefCell = 0;
    pRefValue = 0.0;

    setRefCell
    (
        p(),
        p_rgh(),
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh().needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue - getRefCellValue(p(), pRefCell)
        );
        p_rgh() = p() - rho()*gh();
    }

    mesh().schemes().setFluxRequired(p_rgh().name());
    mesh().schemes().setFluxRequired(alpha1().name());

    // MULES compressed flux is registered in case scalarTransport FO needs it.
    autoPtr<surfaceScalarField> alphaPhiUnPtr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "alphaPhiUn",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", phi().dimensions(), 0.0)
        )
    );
    alphaPhiUnPtr.ptr()->store();

    return;
}

void FoamVOF::createTimeControls()
{

    adjustTimeStep =
     runTime().controlDict().lookupOrDefault("adjustTimeStep", false);

    //Table type entry for maxCo:
    if (runTime().controlDict().found("maxCo"))
    {
        maxCoDataPtr = Function1<scalar>::New("maxCo", runTime().controlDict());
        maxCo = maxCoDataPtr->value(runTime().timeIndex());
        if (maxCoDataPtr->type() == Function1Types::Constant<scalar>::typeName)
        {
            maxCoDataPtr.clear();
        }
    }
    else
    {
        maxCo =
            runTime().controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    }

    maxDeltaT =
        runTime().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    minDeltaT =
        runTime().controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

    return;
}

void FoamVOF::readTimeControls()
{
    adjustTimeStep =
    runTime().controlDict().lookupOrDefault("adjustTimeStep", false);

    //Table type entry for maxCo:
    if (runTime().controlDict().found("maxCo"))
    {
        maxCoDataPtr = Function1<scalar>::New("maxCo", runTime().controlDict());
        maxCo = maxCoDataPtr->value(runTime().timeIndex());
        if (maxCoDataPtr->type() == Function1Types::Constant<scalar>::typeName)
        {
            maxCoDataPtr.clear();
        }
    }
    else
    {
        maxCo =
            runTime().controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    }

    maxDeltaT =
        runTime().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    minDeltaT =
        runTime().controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

    return;
}

void FoamVOF::updateFields
(
    pimpleControl &pimple,
    incompressible::turbulenceModel &turbulence,
    fv::options& fvOptions
)
{
    return;
}

void FoamVOF::setRDeltaT
(
  pimpleControl &pimple
)
{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxAlphaCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    label nAlphaSpreadIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSpreadIter", 1)
    );

    scalar alphaSpreadDiff
    (
    pimpleDict.lookupOrDefault<scalar>("alphaSpreadDiff", 0.2)
    );

    scalar alphaSpreadMax
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    scalar alphaSpreadMin
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMin", 0.01)
    );

    label nAlphaSweepIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSweepIter", 5)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(rhoPhi()))()()
        /((2*maxCo)*mesh().V()*rho())
    );

    if (maxAlphaCo < maxCo)
    {
    // Further limit the reciprocal time-step
    // in the vicinity of the interface

    volScalarField alpha1Bar(fvc::average(alpha1()));

    rDeltaT.ref() = max
    (
        rDeltaT(),
        pos0(alpha1Bar() - alphaSpreadMin)
       *pos0(alphaSpreadMax - alpha1Bar())
       *fvc::surfaceSum(mag(phi()))()()
       /((2*maxAlphaCo)*mesh().V())
    );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
    << gMin(1/rDeltaT.primitiveField())
    << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if (nAlphaSpreadIter > 0)
    {
        fvc::spread
    (
        rDeltaT,
        alpha1(),
        nAlphaSpreadIter,
        alphaSpreadDiff,
        alphaSpreadMax,
        alphaSpreadMin
    );
    }

    if (nAlphaSweepIter > 0)
    {
        fvc::sweep(rDeltaT, alpha1(), nAlphaSweepIter, alphaSpreadDiff);
    }

    Info<< "Smoothed flow time scale min/max = "
    << gMin(1/rDeltaT.primitiveField())
    << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
        && runTime().timeIndex() > runTime().startTimeIndex() + 1
    )
    {
    rDeltaT = max
    (
        rDeltaT,
        (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
    );

    Info<< "Damped flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
    return;
}


void FoamVOF::alphaCourantNoIF()
{

    maxAlphaCo=readScalar(runTime().controlDict().lookup("maxAlphaCo"));
    alphaCoNum = 0.0;
    meanAlphaCoNum = 0.0;

    if (mesh().nInternalFaces())
    {
        scalarField sumPhi
        (
            mixture->nearInterfaceCo()().primitiveField()
            *fvc::surfaceSum(mag(phi()))().primitiveField()
        );

        alphaCoNum = 0.5*gMax(sumPhi/mesh().V().field())*runTime().deltaTValue();

        meanAlphaCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh().V().field()))*runTime().deltaTValue();
    }

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

  return;
}

void FoamVOF::CourantNo()
{
    CoNum = 0.0;
    meanCoNum = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi()))().primitiveField()
    );

    CoNum = 0.5*gMax(sumPhi/mesh().V().field())*runTime().deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh().V().field()))*runTime().deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
    return;
}


void FoamVOF::setInitialDeltaT()
{
  if (adjustTimeStep)
  {
    if ((runTime().timeIndex() == 0) && (CoNum > SMALL))
    {
        const_cast<Time&>(runTime()).setDeltaT
        (
            min
            (
                maxCo*runTime().deltaTValue()/CoNum,
                min(runTime().deltaTValue(), maxDeltaT)
            )
        );
    }

    //exit Foam if deltaT is below user secified minDeltaT
    if (runTime().deltaTValue() < minDeltaT)
    {
        FatalErrorInFunction
            << "deltaT = " <<  runTime().deltaTValue() << nl
            << "minDeltaT = " << minDeltaT <<nl
            << "The computed initial deltaT is lower than minDeltaT"
            << exit(FatalError);
    }

  }
  return;
}

void FoamVOF::setDeltaT()
{

   if (adjustTimeStep)
   {
        scalar maxDeltaTFact =
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        const_cast<Time&>(runTime()).setDeltaT
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

  return;
}

void FoamVOF::correctPhi(pimpleControl &pimple)
{
    CorrectPhi
    (
        U(),
        phi(),
        p_rgh(),
        dimensionedScalar("rAUf", dimTime/rho().dimensions(), 1),
        geometricZeroField(),
        pimple
    );

    continuityErrs();
    return;
}


void FoamVOF::continuityErrs()
{
    volScalarField contErr(fvc::div(phi()));
    scalar meanContErrMag = runTime().deltaTValue()*
    mag(contErr)().weightedAverage(mesh().V()).value();

    scalar globalContErr = runTime().deltaTValue()*
    contErr.weightedAverage(mesh().V()).value();
    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : mean magnitude = " << meanContErrMag
    << ", global = " << globalContErr
    << ", cumulative = " << cumulativeContErr
    << endl;
}


void FoamVOF::alphaEqnSubCycle()
{
    if (cavitation())
    {
        return alphaEqn(alpha1());
    }

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime().deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("0", rhoPhi().dimensions(), 0)
      );

     tmp<volScalarField> trSubDeltaT;

     if (LTS)
     {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh(), nAlphaSubCycles);
     }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1(), nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        alphaEqn();
        rhoPhiSum += (runTime().deltaT()/totalDeltaT)*rhoPhi();
    }

        rhoPhi() = rhoPhiSum;
    }
    else
    {
        alphaEqn();
    }

    rho().forceAssign(alpha1()*rho1() + alpha2()*rho2());
    return;
}


void FoamVOF::alphaControls()
{
    const dictionary& alphaControls = mesh().solution().solverDict(alpha1().name());
    nAlphaSubCycles=1;
    label nAlphaSubCycles0(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    if (nAlphaSubCycles0>1)
    {
        nAlphaSubCycles=nAlphaSubCycles0;
    }
    nAlphaCorr=1;
    label nAlphaCorr0(readLabel(alphaControls.lookup("nAlphaCorr")));
    if (nAlphaCorr0>=0)
    {
        nAlphaCorr=nAlphaCorr0;
    }

    //bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));
    MULESCorr=alphaControls.lookupOrDefault<Switch>("MULESCorr", false);

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed


    alphaApplyPrevCorr=alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false);

// Isotropic compression coefficient
    icAlpha=alphaControls.lookupOrDefault<scalar>("icAlpha", 0);
    return;
}


void FoamVOF::UEqnSolve
(
    fvVectorMatrix& UEqn,
    incompressible::turbulenceModel& turbulence,
    fv::options &fvOptions,
    immiscibleIncompressibleTwoPhaseMixture &mixture,
    pimpleControl &pimple
)
{
    UEqn = fvVectorMatrix
    (
        fvm::ddt(rho(), U())
        + fvm::div(rhoPhi(), U())
        + fvOptions.MRFDDt(rho(), U())
        + turbulence.divDevRhoReff(rho(), U())
    );

    if (cavitation() != 0)
    {
        UEqn -= fvm::Sp(fvc::ddt(rho()) + fvc::div(rhoPhi()), U());
    }

    UEqn -= fvOptions(rho(), U());

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                    - ghf()*fvc::snGrad(rho())
                    - fvc::snGrad(p_rgh())
                ) * mesh().magSf()
            )
        );

        fvOptions.correct(U());
    }
    else
    {
        UEqn.boundaryManipulate(U().boundaryFieldRef());
    }
    return;
}

surfaceScalarField& FoamVOF::ghf()
{
    return mesh_.lookupObjectRef<surfaceScalarField>("ghf");
}


FoamVOF::~FoamVOF()
{
    for (label i=0;i<models.size();++i)
    {
       delete models[i];
    }
}

void FoamVOF::updateAcceleration(autoPtr<Function1<vector>> frameAcceleration)
{

   if (frameAcceleration.valid())
   {
        dimensionedScalar ghRef("ghRef", g().dimensions()*dimLength,ghRefv);
        dimensionedVector fAcc
        (
            "fAcc",
            g().dimensions(),
            frameAcceleration->value(runTime().timeOutputValue())
        );
        gh() = ((g() - fAcc) & mesh().C()) - ghRef;
        ghf() = ((g() - fAcc) & mesh().Cf()) - ghRef;
   }
    return;
}

void FoamVOF::pEqn
(
    fvVectorMatrix& UEqn,
    fv::options &fvOptions,
    immiscibleIncompressibleTwoPhaseMixture &mixture,
    pimpleControl &pimple
)
{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p_rgh()));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
        + fvc::interpolate(rho()*rAU)*fvc::ddtCorr(U(), phi())
    );
    fvOptions.makeRelative(phiHbyA);

    if (mesh().changing())
    {
        if (p_rgh().needReference())
        {
            fvc::makeRelative(phiHbyA, U());
            adjustPhi(phiHbyA, U(), p_rgh());
            fvc::makeAbsolute(phiHbyA, U());
        }
    }
    else
    {
        adjustPhi(phiHbyA, U(), p_rgh());
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
            - ghf()*fvc::snGrad(rho())
        )*rAUf*mesh().magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh(), U(), phiHbyA, rAUf, fvOptions);

    //volume source
    volScalarField su(p_rgh());
    su.primitiveFieldRef()=0.0;

    volScalarField sp(p_rgh());
    sp.primitiveFieldRef()=0.0;

    sp.dimensions().reset(dimensionSet(0,0,-1,0,0,0,0)/dimPressure);
    su.dimensions().reset(dimensionSet(0,0,-1,0,0,0,0));

    if (cavitation())
    {
        forAll(su,i)
        {
            su.primitiveFieldRef()[i]=vdotP().primitiveField()[i]*prop().psat();
            sp.primitiveFieldRef()[i]=vdotP().primitiveField()[i];
        }
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh())
            -fvm::Sp(sp,p_rgh())+fvc::Su(su, p_rgh())
            == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh(), pRefCell));

        p_rghEqn.solve
            (
                mesh().solution().
                    solver(p_rgh().select(pimple.finalInnerIter()))
            );

        if (pimple.finalNonOrthogonalIter())
        {
            phi() = phiHbyA - p_rghEqn.flux();

            p_rgh().relax();

            U() = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U().correctBoundaryConditions();
            fvOptions.correct(U());
        }
   }

    continuityErrs();

    p().forceAssign(p_rgh() + rho()*gh());

    if (p_rgh().needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue - getRefCellValue(p(), pRefCell)
        );
        p_rgh() = p() - rho()*gh();
    }
   return;
}

bool FoamVOF::hasVolSource()
{
    label summdl=0;
    for (label i=0;i<models.size();i++)
    {
        summdl+=models[i]->massTransfer();
    }
    if (summdl>=1)
    {
        return true;
    }
    return false;
}

bool FoamVOF::hasMomentumSource()
{
    label summdl=0;
    for (label i=0;i<models.size();i++)
    {
      summdl+=models[i]->momentumTransfer();
    }
    if (summdl>=1)
    {
        return true;
    }
    return false;
}

bool FoamVOF::hasHeatSource()
{
    label summdl=0;
    for (label i=0;i<models.size();i++)
    {
        summdl+=models[i]->heatTransfer();
    }
    if (summdl>=1)
    {
        return true;
    }
    return false;
}

bool FoamVOF::hasSpecSource()
{
    label summdl=0;
    for (label i=0;i<models.size();i++)
    {
        summdl+=models[i]->specTransfer();
    }
    if (summdl>=1)
    {
        return true;
    }
    return false;
}

void FoamVOF::updateSpec()
{
    for (label i=0;i<models.size();i++)
    {
        models[i]->updateSpec();
    }
    return;
}

void FoamVOF::specSource
(
    volScalarField &Su,
    volScalarField &Sp
)
{
    for (label i=0;i<models.size();i++)
    {
        models[i]->specSource(alpha1().name(),Su,Sp);
    }
    return;
}

void FoamVOF::volSource
(
    volScalarField &Su,
    volScalarField &Sp
)
{
    for (label i=0;i<models.size();i++)
    {
        models[i]->volSource(Su,Sp);
    }
    return;
}

void FoamVOF::alphaEqn(volScalarField &alpha)
{
   // word alphaScheme("div(phi,alpha)");
   // word alpharScheme("div(phirb,alpha)");
    Info<<"Species source term exists..."<<endl;
    volScalarField su(alpha1());
    su=0;
    su.dimensions().reset(dimensionSet(0,0,-1,0,0,0,0));

    volScalarField sp(alpha1());
    sp=0;
    sp.dimensions().reset(dimensionSet(0,0,-1,0,0,0,0));

    specSource(su,sp);

    fvScalarMatrix alphaEqn
    (
       fvm::ddt(alpha)
       + fvm::div(phi(), alpha)
       +fvm::Sp(sp,alpha)
       - fvc::Su(su,alpha)
    );

    alphaEqn.relax();
    alphaEqn.solve();

    alpha2() = 1.0 - alpha;
    updateRho();
    updateRhoPhi();

    return;
}

void FoamVOF::alphaEqn()
{

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
    tmp<fv::ddtScheme<scalar>> tddtAlpha
    (
        fv::ddtScheme<scalar>::New
        (
            mesh(),
            mesh().schemes().ddtScheme("ddt(alpha)")
        )
    );
    const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

    if
    (
        isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
        || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
    )
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
            << "Sub-cycling is not supported "
               "with the CrankNicolson ddt scheme"
            << exit(FatalError);
        }

        if
        (
            alphaRestart
            || mesh().time().timeIndex() > mesh().time().startTimeIndex() + 1
        )
        {
            ocCoeff =
            refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
           .ocCoeff();
        }
    }
    else
    {
        FatalErrorInFunction
        << "Only Euler and CrankNicolson ddt schemes are supported"
        << exit(FatalError);
    }
}

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture->cAlpha()*mag(phi()/mesh().magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture->cAlpha()*icAlpha)*fvc::interpolate(mag(U()));
    }

    surfaceScalarField::Boundary& phicBf =
    phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp.forceAssign(0);
        }
    }

    tmp<surfaceScalarField> phiCN(phi());

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi() + (1.0 - cnCoeff)*phi().oldTime();
    }

    if (MULESCorr)
    {

        zeroField Su;
        zeroField Sp;
        zeroField divU;

        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
                ? fv::localEulerDdtScheme<scalar>(mesh()).fvmDdt(alpha1())
                : fv::EulerDdtScheme<scalar>(mesh()).fvmDdt(alpha1())
            )
            + fv::gaussConvectionScheme<scalar>
            (
                mesh(),
                phiCN,
                upwind<scalar>(mesh(), phiCN)
            ).fvmDiv(phiCN, alpha1())
            ==
            Su + fvm::Sp(Sp + divU, alpha1())
        );

    alpha1Eqn.solve();

    Info<< "Phase-1 volume fraction = "
        << alpha1().weightedAverage(mesh().Vsc()).value()
        << "  Min(" << alpha1().name() << ") = " << min(alpha1()).value()
        << "  Max(" << alpha1().name() << ") = " << max(alpha1()).value()
        << endl;

    tmp<surfaceScalarField> talphaPhiUD(alpha1Eqn.flux());
    alphaPhi() = talphaPhiUD();

    if (alphaApplyPrevCorr && talphaPhiCorr0.valid())
    {
        Info<< "Applying the previous iteration compression flux" << endl;
        MULES::correct(alpha1(), alphaPhi(), talphaPhiCorr0.ref(), 1, 0);

        alphaPhi() += talphaPhiCorr0();
    }

        // Cache the upwind-flux
        talphaPhiCorr0 = talphaPhiUD;

        alpha2() = 1.0 - alpha1();

        mixture->correct();
    }


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        zeroField Su;
        zeroField Sp;
        zeroField divU;

        surfaceScalarField phir(phic*mixture->nHatf());

    tmp<surfaceScalarField> talphaPhiUn
    (
        fvc::flux
        (
            phiCN(),
            cnCoeff*alpha1() + (1.0 - cnCoeff)*alpha1().oldTime(),
            alphaScheme
        )
        + fvc::flux
        (
           -fvc::flux(-phir, alpha2(), alpharScheme),
            alpha1(),
            alpharScheme
        )
    );

    if (MULESCorr)
    {
        tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi());
        volScalarField alpha10("alpha10", alpha1());

        MULES::correct
        (
            geometricOneField(),
            alpha1(),
            talphaPhiUn(),
            talphaPhiCorr.ref(),
            Sp,
            (-Sp*alpha1())(),
            1,
            0
        );

        // Under-relax the correction for all but the 1st corrector
        if (aCorr == 0)
        {
            alphaPhi() += talphaPhiCorr();
        }
        else
        {
            alpha1() = 0.5*alpha1() + 0.5*alpha10;
            alphaPhi() += 0.5*talphaPhiCorr();
        }
    }
    else
    {
        alphaPhi() = talphaPhiUn;

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1(),
            phiCN,
            alphaPhi(),
            Sp,
            (Su + divU*min(alpha1(), scalar(1)))(),
            1,
            0
        );
    }

    alpha2() = 1.0 - alpha1();

    mixture->correct();
  }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhiCorr0 = alphaPhi() - talphaPhiCorr0;
        talphaPhiCorr0.ref().rename("alphaPhiCorr0");
    }
    else
    {
        talphaPhiCorr0.clear();
    }

    if
    (
        word(mesh().schemes().ddtScheme("ddt(rho,U)"))
        == fv::EulerDdtScheme<vector>::typeName
    )
    {
        const dimensionedScalar& rho1f(rho1());
        const dimensionedScalar& rho2f(rho2());

        rhoPhi() = alphaPhi()*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi() = (alphaPhi() - (1.0 - cnCoeff)*alphaPhi().oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        const dimensionedScalar& rho1f(rho1());
        const dimensionedScalar& rho2f(rho2());
        rhoPhi() = alphaPhi()*(rho1f - rho2f) + phi()*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
    << alpha1().weightedAverage(mesh().Vsc()).value()
    << "  Min(" << alpha1().name() << ") = " << min(alpha1()).value()
    << "  Max(" << alpha1().name() << ") = " << max(alpha1()).value()
    << endl;


  return;
}

//material properties update

void FoamVOF::updateRho()
{

  rho()=rho1()*alpha1()+rho2()*alpha2();

  return;
}

void FoamVOF::updateRhoPhi()
{
    surfaceScalarField alphaf(fvc::interpolate(alpha1()));
    rhoPhi()=(rho1()-rho2())*alphaf*phi()+rho2()*phi();
   // Info<<"rho1,rho2:"<<rho1()<<","<<rho2()<<endl;
}

void FoamVOF::updateCp()
{
  if (solveT()==0) return;
  dimensionedScalar cp_1("cp", dimensionSet(0,2,-2,-1,0,0,0),prop().cp1());
  dimensionedScalar cp_2("cp",dimensionSet(0,2,-2,-1,0,0,0) ,prop().cp2());

  dimensionedScalar rho_1("rho", dimensionSet(1,-3,0,0,0,0,0),prop().rho1());
  dimensionedScalar rho_2("rho", dimensionSet(1,-3,0,0,0,0,0),prop().rho2());
  cp()=((cp_1*rho_1)*alpha1()+(cp_2*rho_2)*alpha2())/rho();

  return;
}

void FoamVOF::updateK()
{
  if (solveT()==0) return;
  dimensionedScalar k_1("K", dimensionSet(1,1,-3,-1,0,0,0),prop().k1());
  dimensionedScalar k_2("K",dimensionSet(1,1,-3,-1,0,0,0) ,prop().k2());

  K()=k_1*alpha1()+k_2*alpha2();

  return;
}

void FoamVOF::TEqnSolve
(
   incompressible::turbulenceModel& turbulence
)
{
    scalar prv=1.0/prop().Pr();

    fvScalarMatrix TEqn
    (
      fvm::ddt(T())
      + fvm::div(phi(), T())

      - 1.0/(rho()*cp())*
      (
    fvm::laplacian
    (
     K()+prv*turbulence.nut()*rho()*cp(),T()
    )
      )
    );
    TEqn -= 1.0/(rho()*cp())*fvc::laplacian(prv*turbulence.nut()*rho()*T(),cp());
    TEqn.relax();
    TEqn.solve();

    return;
}

void FoamVOF::dynamicMeshUpdate
(
    pimpleControl &pimple
)
{
    scalar timeBeforeMeshUpdate = runTime().elapsedCpuTime();

    const_cast<dynamicFvMesh&>(mesh()).update();

    if (mesh().changing())
    {
        Info<< "Execution time for mesh.update() = "
            << runTime().elapsedCpuTime() - timeBeforeMeshUpdate
            << " s" << endl;

        // Do not apply previous time-step mesh compression flux
        // if the mesh topology changed
        if (mesh().topoChanging())
        {
            talphaPhiCorr0.clear();
        }

        dimensionedScalar ghRef("ghRef", g().dimensions()*dimLength,ghRefv);
        gh() = (g() & mesh().C()) - ghRef;
        ghf() = (g() & mesh().Cf()) - ghRef;
    }

    if (mesh().changing() || mesh().topoChanging())
    {
        // Calculate absolute flux from the mapped surface velocity
        // Note: temporary fix until mapped Uf is assessed
        Uf() = fvc::interpolate(U());

        // Calculate absolute flux from the mapped surface velocity
        phi() = mesh().Sf() & Uf();

        correctPhi(pimple);

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi(), U());

        mixture->correct();
    }

    //if (mesh.changing() && checkMeshCourantNo)
    //{
    //    #include "include/meshCourantNo.H"
    //}
}
