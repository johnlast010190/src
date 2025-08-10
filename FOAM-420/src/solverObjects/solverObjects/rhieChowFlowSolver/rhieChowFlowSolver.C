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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "rhieChowFlowSolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/CourantNumber/CourantNumber.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "cfdTools/general/adjustSplitBCs/adjustSplitBCs.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(rhieChowFlowSolver, 0);
}
}

makeFvSolverOption(rhieChowFlowSolver);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::rhieChowFlowSolver::rhieChowFlowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    flowSolver(name, obr, dict),
    thermoPtr_(nullptr),
    ddtPhiCorr_(true),
    printContErr_(true),
    solnControlPtr_(nullptr),
    useGradP_(false),
    transient_(false),
    buoyant_(false),
    distinctBuoyancy_(false)
{
    const Time& runTime = mesh().time();

    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    transient_ = !isA<simpleControl>(*solnControlPtr_);

    if (!obr_.foundObject<volScalarField>("p"))
    {
        Info<< "Creating initial p field\n" << endl;
        p_ = new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );
        obr_.store(p_);
    }
    else
    {
        Info<< "Reading field p\n" << endl;
        p_ = obr.lookupObjectRefPtr<volScalarField>("p");
    }
    volScalarField& p(*p_);

    Info<< "Reading field U\n" << endl;
    U_.set
    (
        new volVectorField
        (
            IOobject
            (
                IOobject::groupName("U", phaseName_),
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<rhoThermo>(multiphaseThermo::lookupOrCreate(obr, phaseName_));

    // Initialize info from material models
    buoyant_ = thermo().buoyant();
    distinctBuoyancy_ = thermo().distinctBuoyancy();

    Info<< "Creating field rho\n" << endl;
    rho_.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("rho", phaseName_),
                runTime.timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo().rho()
        )
    );

    if (buoyant_)
    {
        Info<< "Creating buoyant rho\n" << endl;
        bRho_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("buoyantRho", phaseName_),
                    runTime.timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_()
            )
        );
        if (distinctBuoyancy_)
        {
            bRho_() = thermo().buoyantRho();
        }

        initializeGravityHref();
    }

    updateSettingsFromDict();

    setRefCell
    (
        p,
        solnControlPtr_->dict(),
        pRefCell_,
        pRefValue_ //From solution dict
    );
    if (!thermo().incompressible())
    {
        compressibility_ =
            new dimensionedScalar(fvc::domainIntegrate(thermo().psi()));
    }
    if (p.needReference() && thermo().incompressible())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue_ - getRefCellValue(p, pRefCell_)
        );
    }

    if (verbose_)
    {
        Info<< "p min/max: "
            << min(p).value() << " "
            << max(p).value() << endl;
        Info<< "rho min/max: "
            << min(rho_()).value() << " "
            << max(rho_()).value() << endl;

        if (distinctBuoyancy_)
        {
            Info<< "Buoyant rho min/max: "
                << min(bRho_()).value() << " "
                << max(bRho_()).value() << endl;
        }
    }

    mesh().schemes().setFluxRequired(p.name());

    initialMass_.set(new dimensionedScalar(fvc::domainIntegrate(rho_())));

    rhoMin_.set
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMin",
                solnControlPtr_->dict(),
                dimDensity,
                0
            )
        )
    );

    rhoMax_.set
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMax",
                solnControlPtr_->dict(),
                dimDensity,
                GREAT
            )
        )
    );

    // Read or create face flux
    IOobject phiIOObj
    (
        IOobject::groupName("phi", phaseName_),
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phiIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating face flux field phi" << nl << endl;
        recreatePhi_ = true;
    }
    else
    {
        Info<< "Reading face flux field phi" << nl << endl;
        recreatePhi_ = false;
    }
    phi_.set
    (
        new surfaceScalarField
        (
            phiIOObj,
            linearInterpolate(rho_()*U_()) & mesh().Sf()
        )
    );

    if (solnControlPtr_->transonic())
    {
        phiv_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phiv", phaseName_),
                    runTime.timeName(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(U_()) & mesh().Sf()
            )
        );
        // Make flux consistent with density interpolation
        phi_() = rhoInterpolation()*phiv_();
    }
    if (!isStatic())
    {
        Info<< "Reading/calculating field rhoUf\n" << endl;
        rhoUf_.set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName("rhoUf", phaseName_),
                    runTime.timeName(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho_()*U_())
            )
        );
        fvc::makeRelative(phi_(), rho_(), U_());
    }

    pressureControl_.set
    (
        new pressureControl
        (
            p,
            rho_,
            solnControlPtr_->dict(),
            thermo().pRefValue()
        )
    );

    cumulativeContErr_ = 0;


    Info<< "Creating turbulence model\n" << endl;
    turbulence_=
    (
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            phi_,
            thermo()
        ).ptr()
    );
    turbulence_->store();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rhieChowFlowSolver::updateSettingsFromDict()
{
    ddtPhiCorr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "ddtPhiCorr",
        true
    );

    printContErr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "printContinuityErrors",
        true
    );

    useGradP_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "useGradP",
        !buoyant_
    );
}

void Foam::fv::rhieChowFlowSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    word UPredName = U_->name() + "Predictor";
    solveNames = {UPredName, p_->name(), U_->name()};

    // Dependencies
    optionalDependencies.insert(UPredName, {"fvMesh"});
    requiredDependencies.insert(p_->name(), {UPredName});
    requiredDependencies.insert(U_->name(), {p_->name()});

    // Correctors

    // Add all solves to outer corrector
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);

    // PISO (inner) correctors
    correctorMembers.insert("PISOCorrector", {p_->name(), U_->name()});

    // Non-orthogonal correctors
    correctorMembers.insert
    (
        "nonOrthogonalCorrector:" + p_->name(), {p_->name()}
    );

    optionalDependencies.insert("materialProperties", {p_->name()});
}


bool Foam::fv::rhieChowFlowSolver::initialise()
{
    volScalarField& p = *p_;

    if (recreatePhi_ && !transient_)
    {
        tmp<fvScalarMatrix> fvOptSrc
        (
            this->fvOptions()
            (
                thermo().psi(), p, rho_().name()
            )
        );

        CorrectPhi
        (
            U_(),
            phi_(),
            p,
            rho_(),
            thermo().psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            (fvOptSrc & p)(),
            *solnControlPtr_,
            p.name()
        );
        // updateCoeffs was called in matrix construction but evaluate was
        // not called - resetUpdate reverts to updatable state
        forAll(p.boundaryField(), patchi)
        {
            p.boundaryFieldRef()[patchi].resetUpdate();
        }

        Info<< endl;
    }

    if (!exists(phi_().objectPath()))
    {
        tmp<surfaceScalarField> rhof(fvc::interpolate(rho_()));
        this->fvOptions().makeRelative(rhof, phi_());
    }

    if (transient_)
    {
        // Calculate initial Courant number
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());
    }

    return true;
}


scalar Foam::fv::rhieChowFlowSolver::getMaxTimeStep()
{
    // Try region-specific maxCo, otherwise revert to global
    const dictionary& maxCoDict
    (
        solnControlPtr_->dict().found("maxCo")
      ? solnControlPtr_->dict()
      : mesh_.time().controlDict()
    );
    maxCoDataPtr_ = Function1<scalar>::New("maxCo", maxCoDict);

    scalar maxMaxCo = maxCoDataPtr_->value(mesh_.time().timeIndex());
    return maxMaxCo/stabilise(maxCo_, SMALL)*mesh_.time().deltaTValue();
}


const Foam::volScalarField& Foam::fv::rhieChowFlowSolver::buoyantRho() const
{
    return bRho_();
}


bool Foam::fv::rhieChowFlowSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if (correctorName == "PISOCorrector")
    {
        if (transient_)
        {
            label nCorr = refCast<pimpleControl>(*solnControlPtr_).nCorrPISO();
            return (corrector >= nCorr-1);
        }
        else
        {
            return true;
        }
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return true;
    }
}


void Foam::fv::rhieChowFlowSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    // Update every time step
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        updateSettingsFromDict();
    }

    if (correctorName == solverObject::outerCorrectorName)
    {
        if (!isStatic())
        {
            // Store divrhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
            divRhoU_.clear();
            divRhoU_ =
                tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "divrhoU",
                        fvc::div(fvc::absolute(phi_(), rho_(), U_()))
                    )
                );

            // Store momentum to set rhoUf for introduced faces.
            rhoU_.clear();
            rhoU_ =
                tmp<volVectorField>(new volVectorField("rhoU", rho_()*U_()));
        }
    }
    else if (correctorName == "PISOCorrector")
    {
        corr_ = corrector;
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        nonOrthCorr_ = corrector;
    }
}


tmp<fvVectorMatrix>
Foam::fv::rhieChowFlowSolver::assembleVectorMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if
    (
        fieldName == U_->name()+"Predictor"
     && solnControlPtr_->momentumPredictor()
    )
    {
        if (correctorNumber_[solverObject::outerCorrectorName] == 0)
        {
            // Explicit predictor for new value of rho
            // Only use if compressible and running in PISO mode
            if
            (
                transient_
             && !thermo().incompressible()
             && solnControlPtr_->nOuterCorr() <= 1
            )
            {
                solveContinuity();
            }
            rho_().storePrevIter();
        }

        return UPredictorEqn();
    }
    else
    {
        return tmp<fvVectorMatrix>();
    }
}


tmp<fvVectorMatrix> Foam::fv::rhieChowFlowSolver::UPredictorEqn()
{
    assembleUEqnLHS();

    tmp<fvVectorMatrix> tUEqn(new fvVectorMatrix(UEqnLHS_()));
    if (useGradP_)
    {
        tUEqn.ref() += fvc::grad(*p_);
        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce();
            tUEqn.ref() -=
                fvc::reconstruct(faceBuoyancyForce_()*mesh_.magSf());
        }
    }
    else
    {
        tmp<surfaceScalarField> tsnGrad;
        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce(true);
            tsnGrad = -faceBuoyancyForce_();
        }
        else
        {
            tsnGrad = fvc::snGrad(*p_);
        }
        tUEqn.ref() += fvc::reconstruct(tsnGrad*mesh_.magSf());
    }
    return tUEqn;
}


void Foam::fv::rhieChowFlowSolver::assembleUEqnLHS()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();

    fv::options& fvOptions = this->fvOptions();

    if (!transient_ && solnControlPtr_->modifiedMomentumInterp())
    {
        U.storePrevIter();
    }

    // Assemble the Momentum equation
    UEqnLHS_ =
        new fvVectorMatrix
        (
            fvm::div(phi, U)
          + fvOptions.MRFDDt(rho, U)
          + turbulence_->divDevRhoReff(U)
        ==
            fvOptions(rho, U)
        );

    if (transient_)
    {
        if (!thermo().isochoric() || !isStatic())
        {
            UEqnLHS_.ref() += fvm::ddt(rho, U);
        }
        else
        {
            const word schemeName("ddt(" + rho.name() + "," + U.name() + ")");
            UEqnLHS_.ref() += rho*fvm::ddt(U, schemeName);
        }
    }

    UEqnLHS_->relax();

    fvOptions.constrain(UEqnLHS_.ref());
}


tmp<fvScalarMatrix>
Foam::fv::rhieChowFlowSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if (fieldName == p_->name())
    {
        assemblepEqn();
        finalSolve =
            (
                transient_
             && finalIter_[solverObject::outerCorrectorName]
             && finalIter_["PISOCorrector"]
             && finalIter_["nonOrthogonalCorrector:"+p_->name()]
            );
        // Return copy of the tmp so it doesn't get deleted from outside
        return tmp<fvScalarMatrix>(pEqn_);
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }

}


void Foam::fv::rhieChowFlowSolver::assemblepEqn()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();
    volScalarField& p = *p_;
    const volScalarField& psi = thermo().psi();
    fv::options& fvOptions = this->fvOptions();

    if
    (
        !solnControlPtr_->momentumPredictor()
     && corr_ == 0
     && nonOrthCorr_ == 0
    )
    {
        assembleUEqnLHS();
    }

    if (nonOrthCorr_ == 0)
    {
        if (!thermo().isochoric())
        {
            rho = thermo().rho();

            // Relax density here in response to temperature change. Relaxed again
            // in correct() to respond to pressure change
            // Store for consistency in case no relaxation in use
            rho.relax();
        }
        if (buoyant_)
        {
            if (distinctBuoyancy_)
            {
                bRho_() = thermo().buoyantRho();
            }
            else
            {
                bRho_() = rho;
            }
        }

        if (!thermo().incompressible())
        {
            // Thermodynamic density needs to be updated by psi*d(p) after the
            // pressure solution
            p0_ = tmp<volScalarField>(new volScalarField(p));
        }

        tmp<volScalarField> trAU
        (
            new volScalarField(1/UEqnLHS_().A())
        );

        // Fishy treatment inside constranHbyA (needs to be looked at)
        // When used directly HbyA_ = constrainHbyA(trAU()*UEqnLHS_().H(), U, p);
        // The tmp<> seems to get corrupted
        HbyA_ = new volVectorField(constrainHbyA(trAU()*UEqnLHS_().H(), U, p));
        volVectorField& HbyA = HbyA_.ref();
        if (solnControl().consistent())
        {
            tmp<volScalarField> A(1/trAU());
            tmp<volScalarField> AmH1(-UEqnLHS_().H1()+A());
            // Stabilise the case where H1 ~= -A
            rAtU_ = new volScalarField("(1|At(U))", 1/max(AmH1, 0.01*A));
        }
        else
        {
            // Transfer the tmp
            rAtU_ = trAU;
        }
        const volScalarField& rAtU = rAtU_();

        tmp<surfaceScalarField> phivHbyA
        (
            new surfaceScalarField("phivHbyA", fvc::flux(HbyA))
        );
        tmp<surfaceScalarField> rAtUf
        (
            new surfaceScalarField
            (
                "rAUf",
                fvc::interpolate(rAtU, "interpolate((1|A(U)))")
            )
        );

        setOrComputeRhof();

        if (solnControlPtr_->transonic())
        {
            rhorAtUf_ = new surfaceScalarField("rhorAUf", rhof_()*rAtUf());
            phiHbyA_ = new surfaceScalarField("phiHbyA", rhof_()*phivHbyA());
        }
        else
        {
            // Reuse tmps for efficiency
            rhorAtUf_ = new surfaceScalarField("rhorAUf", rhof_()*rAtUf);
            phiHbyA_ = new surfaceScalarField("phiHbyA", rhof_()*phivHbyA);
        }
        const surfaceScalarField& rhorAtUf = rhorAtUf_();
        surfaceScalarField& phiHbyA = phiHbyA_.ref();

        tmp<surfaceScalarField> rhorAUf;
        if (solnControlPtr_->consistent())
        {
            tmp<surfaceScalarField> rAUf =
                fvc::interpolate(trAU(), "interpolate((1|A(U)))");
            rhorAUf = rhof_()*rAUf;
        }
        else
        {
            // Point tmp to existing var; increments ref count
            rhorAUf = tmp<surfaceScalarField>(rhorAtUf_);
        }

        if (transient_ && ddtPhiCorr_)
        {
            if (rhoUf_.valid())
            {
                //- Dynamic mesh
                phiHbyA +=
                    fvOptions.zeroFilter
                    (
                        rhorAUf()*fvc::ddtCorr(rho, U, rhoUf_())
                    );
            }
            else
            {
                //- Static mesh
                phiHbyA +=
                    fvOptions.zeroFilter
                    (
                         rhorAUf()*fvc::ddtCorr(rho, U, phi)
                    );
            }
        }
        fvOptions.makeRelative(rhof_(), phiHbyA);

        // If phi was initialised with correctPhi, this causes a large
        // discrepancy here between phi and U, so we skip it for the first iter
        if
        (
            solnControlPtr_->modifiedMomentumInterp()
         && !transient_
         && !recreatePhi_
        )
        {
            tmp<surfaceScalarField> rho0f =
                fvc::interpolate(rho.prevIter(), "interpolate(rho)");
            fvOptions.makeAbsolute(rho0f(), phi);
            phiHbyA +=
                (1 - mesh_.solution().equationRelaxationFactor(U.name()))
               *rho0f()
               *(
                    fvc::absolute(phi/rho0f(), U)
                  - (fvc::interpolate(U.prevIter()) & mesh_.Sf())
                );
        }
        recreatePhi_ = false;

        if (transient_)
        {
            closedVolume_ = p.needReference();
        }
        else
        {
            closedVolume_ = adjustPhi(phiHbyA, U, p);
        }

        if (solnControl().consistent())
        {
            tmp<surfaceScalarField> tsnGrad;
            if (useGradP_)
            {
                tsnGrad = fvc::snGrad(p)*mesh().magSf();
                tmp<volVectorField> tgrad(fvc::grad(p));
                if (buoyant_)
                {
                    tsnGrad.ref() -= faceBuoyancyForce_()*mesh_.magSf();
                    tgrad.ref() -=
                        fvc::reconstruct(faceBuoyancyForce_*mesh_.magSf());
                }
                phiHbyA += (rhorAtUf-rhorAUf())*tsnGrad();
                HbyA -= (trAU() - rAtU)*tgrad;
            }
            else
            {
                if (buoyant_)
                {
                    tsnGrad = -faceBuoyancyForce_*mesh().magSf();
                }
                else
                {
                    tsnGrad = fvc::snGrad(p)*mesh().magSf();
                }
                phiHbyA += (rhorAtUf-rhorAUf())*tsnGrad();
                HbyA -= (trAU() - rAtU)*fvc::reconstruct(tsnGrad);
            }
        }

        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce();
            phiHbyA += rhorAtUf_()*faceBuoyancyForce_()*mesh_.magSf();
        }

        // Update the fixedFluxPressure BCs to ensure flux consistency
        constrainPressure(p, rho, U, phiHbyA, rhorAtUf_(), fvOptions);

        if (!thermo().incompressible())
        {
            compressibility_() = fvc::domainIntegrate(psi);
        }

        // For correct behaviour this term has to be after constrain pressure
        // It compensates the fvc::ddt(rho) that produce additional flux (mesh)
        // Even for constant rho
        if (!isStatic())
        {
            fvc::makeRelative(phiHbyA, rho, U);
        }

        pDDtEqn_ =
            new fvScalarMatrix
            (
              - fvOptions(psi, p, rho.name())
              + fvc::div(phiHbyA)
            );
        if (!thermo().isochoric() || !isStatic())
        {
            pDDtEqn_.ref() += fvc::ddt(rho);
        }
        if (!thermo().incompressible())
        {
            pDDtEqn_.ref() += psi*correction(fvm::ddt(p));
        }

        if (solnControlPtr_->transonic())
        {
            phiv_() =
                phivHbyA-rAtUf*fvc::snGrad(p)*mesh_.magSf();
            gaussConvectionScheme<scalar> convScheme
            (
                mesh_,
                phiv_(),
                tmp<surfaceInterpolationScheme<scalar>>
                (
                    new upwind<scalar>(mesh_, phiv_())
                )
            );
            pDDtEqn_.ref() += correction(convScheme.fvmDiv(phiv_(), psi, p));
        }
    }

    // Solve pressure
    pEqn_ =
        new fvScalarMatrix
        (
            pDDtEqn_()
          - fvm::laplacian(rhorAtUf_(), p)
        );

    if (solnControlPtr_->transonic())
    {
        pEqn_->relax();
    }

    // Do not set reference in transient, compressible case as there is a dp/dt
    // term
    if (p.needReference() && (!transient_ || thermo().incompressible()))
    {
        // In the steady, compressible case do not reset to ref. value as the
        // pressure is later modified to respect global mass conservation
        pEqn_->setReference
        (
            pRefCell_,
            thermo().incompressible()
          ? pRefValue_
          : getRefCellValue(p, pRefCell_)
        );
    }
}


void Foam::fv::rhieChowFlowSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();
    if (solveName == U_->name()+"Predictor")
    {
        fvOptions.correct(U_());

        if (verbose_)
        {
            Info<< "U predictor: mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }
    }
    else if (solveName == p_->name())
    {
        if (verbose_)
        {
            Info<< "p min/max: "
                << min(*p_).value() << " "
                << max(*p_).value() << endl;
        }

        if (finalIter_["nonOrthogonalCorrector:"+p_->name()])
        {
            // Final non-orthogonal corrector

            volScalarField& p = *p_;
            volScalarField& rho = rho_();

            // Calculate the conservative flux
            phi_() = phiHbyA_ + pEqn_->flux();

            adjustSplitBCs(phi_(), U_(), p);
            adjustSplitBCs(phi_(), U_(), U_());
            p.relax();
            pressureControl_->limit(p);

            if (printContErr_)
            {
                if (!thermo().incompressible())
                {
                    if (transient_)
                    {
                        // Here this is only used for printing continuity errors
                        // as rho is overwritten below
                        solveContinuity();
                    }
                    printCompContinuityErrors(regionName);
                }
                else
                {
                    printIncoContinuityErrors(regionName);
                }
            }

            if (closedVolume_)
            {
                if (thermo().incompressible())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pRefValue_ - getRefCellValue(p, pRefCell_)
                    );
                }
                else
                {
                    // Ensure that density correction conserves
                    // initial mass
                    p +=
                        (
                            initialMass_()
                          - fvc::domainIntegrate
                            (
                                thermo().rho() + thermo().psi()*(p-p0_())
                            )
                        )/compressibility_();
                }
            }

            if (!thermo().incompressible())
            {
                thermoPtr_->correctRho(thermo().psi()*(p - p0_()));
            }

            if (!thermo().isochoric())
            {
                rho = thermo().rho();
                if (!transient_)
                {
                    rho = max(rho, rhoMin_());
                    rho = min(rho, rhoMax_());
                }
                rho.relax();
            }
            else
            {
                // Keep time index up to date
                rho.timeIndex() = mesh().time().timeIndex();
            }
        }
    }
    else if (solveName == U_->name())
    {
        // U corrector

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        volVectorField& U = U_();
        volScalarField& p = *p_;
        volScalarField& rho = rho_();
        surfaceScalarField& phi = phi_();
        if (useGradP_)
        {
            U = HbyA_ - rAtU_()*fvc::grad(p);
            if (buoyant_)
            {
                U += rAtU_()*fvc::reconstruct(faceBuoyancyForce_*mesh_.magSf());
            }
        }
        else
        {
            //tmp<surfaceScalarField> pEqnFlux(pEqn_->flux()/rhorAtUf_());
            // At present we can't use flux here because pressure relaxation
            // may affect pressure boundaries, yet this is frozen into the
            // boundaryCoeffs so will not reflect in flux
            tmp<surfaceScalarField> pEqnFlux = -fvc::snGrad(p);
            if (buoyant_)
            {
                pEqnFlux.ref() += faceBuoyancyForce_();
            }
            pEqnFlux.ref() *= mesh_.magSf();
            U = HbyA_ + rAtU_()*fvc::reconstruct(pEqnFlux);
        }

        pEqn_.clear();

        U.correctBoundaryConditions();
        fvOptions.correct(U);

        if (verbose_)
        {
            Info<< "mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }

        if (!isStatic())
        {
            rhoUf_() = fvc::interpolate(rho*U);
            surfaceVectorField n(mesh_.Sf()/mesh_.magSf());
            rhoUf_() +=
                n
                *(
                    fvc::absolute(phi, rho, U)/mesh_.magSf()
                    - (n & rhoUf_())
                );
        }

        if (finalIter_["PISOCorrector"] && transient_)
        {
            if (!thermo().isochoric())
            {
                // To undo the relaxation?
                rho = thermo().rho();
            }
        }
    }
}


void Foam::fv::rhieChowFlowSolver::solveContinuity()
{
    fv::options& fvOptions = this->fvOptions();
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_())
      + fvc::div(phi_())
      ==
        fvOptions(rho_())
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho_());
}


void Foam::fv::rhieChowFlowSolver::printIncoContinuityErrors
(
    const word& regionName
)
{
    volScalarField contErr(fvc::div(phi_()));

    scalar sumLocalContErr = mesh_.time().deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    scalar globalContErr = mesh_.time().deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value();
    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::rhieChowFlowSolver::printCompContinuityErrors
(
    const word& regionName
)
{
    const volScalarField& rho = rho_();

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr, globalContErr;

    if (transient_)
    {
        sumLocalContErr =
        (
            fvc::domainIntegrate(mag(rho - thermo().rho()))/totalMass
        ).value();

        globalContErr =
        (
            fvc::domainIntegrate(rho - thermo().rho())/totalMass
        ).value();

        cumulativeContErr_ += globalContErr;
    }
    else
    {
        dimensionedScalar totalMass = fvc::domainIntegrate(rho);

        fv::options& fvOptions = this->fvOptions();
        volScalarField rhoSource ( fvOptions(rho_()) & rho );
        sumLocalContErr =
        (
            fvc::domainIntegrate
            (
                mag(fvc::div(phi_())-rhoSource)
            )/totalMass
        ).value();

        globalContErr =
        (
            fvc::domainIntegrate(fvc::div(phi_())-rhoSource)/totalMass
        ).value();

        cumulativeContErr_ += globalContErr;
    }

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::rhieChowFlowSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if
    (
        transient_
     && correctorName == solverObject::outerCorrectorName
     && finalIter
    )
    {
        Info<< endl;
        // Calculate Courant for next time step
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());
    }
}


void Foam::fv::rhieChowFlowSolver::correctPhi(bool defaultCorrectPhi)
{
    if
    (
        solnControlPtr_->dict().lookupOrDefault<Switch>
        (
            "correctPhi", defaultCorrectPhi
        )
    )
    {
        // Calculate absolute flux from the mapped surface velocity
        phi_() = mesh_.Sf() & rhoUf_();

        CorrectPhi
        (
            U_(),
            phi_(),
            *p_,
            rho_(),
            thermo().psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            divRhoU_(),
            *solnControlPtr_
        );

        // Make the fluxes relative to the mesh-motion
        fvc::makeRelative(phi_(), rho_(), U_());
    }
}


bool Foam::fv::rhieChowFlowSolver::movePoints()
{
    calculateghFields(buoyant_);
    return true;
}


void Foam::fv::rhieChowFlowSolver::updateMesh(const mapPolyMesh& mpm)
{
    correctPhi(true);
    calculateghFields(buoyant_);
}


// ************************************************************************* //
