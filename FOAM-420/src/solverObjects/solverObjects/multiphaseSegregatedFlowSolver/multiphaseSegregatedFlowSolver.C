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
    (c) 2019-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseSegregatedFlowSolver.H"
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
#include "rhoThermo/rhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseSegregatedFlowSolver, 0);
}
}

makeFvSolverOption(multiphaseSegregatedFlowSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseSegregatedFlowSolver::multiphaseSegregatedFlowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    flowSolver(name, obr, dict),
    thermoPtr_(nullptr),
    printContErr_(true),
    solnControlPtr_(nullptr),
    buoyant_(false)
{
    const Time& runTime = mesh().time();

    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    printContErr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "printContinuityErrors",
        true
    );

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
                "U",
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ = &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr));

    buoyant_ = thermoPtr_->buoyant();
    if (thermoPtr_->distinctBuoyancy())
    {
        FatalErrorInFunction
            << "This solver does not support distinct buoyancy." << nl << endl;
    }

    Info<< "Creating field rho\n" << endl;
    rho_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoPtr_->rho()
        )
    );

    if (buoyant_)
    {
        initializeGravityHref();
    }

    setRefCell
    (
        p,
        solnControlPtr_->dict(),
        pRefCell_,
        pRefValue_ //From solution dict
    );
    if (!thermoPtr_->incompressible())
    {
        compressibility_ =
            new dimensionedScalar(fvc::domainIntegrate(thermoPtr_->psi()));
    }
    if (p.needReference() && thermoPtr_->incompressible())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue_ - getRefCellValue(p, pRefCell_)
        );
    }

    mesh().schemes().setFluxRequired(p.name());

    initialMass_.set(new dimensionedScalar(fvc::domainIntegrate(rho_())));

    // Read or create volumetric face flux
    IOobject phivIOObj
    (
        "phiv",
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phivIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating volumetric face flux field phiv" << nl << endl;
        recreatePhi_ = true;
    }
    else
    {
        Info<< "Reading volumetric face flux field phiv" << nl << endl;
        recreatePhi_ = false;
    }
    phiv_.set
    (
        new surfaceScalarField
        (
            phivIOObj,
            linearInterpolate(U_()) & mesh().Sf()
        )
    );

    // Read or create face flux
    IOobject phiIOObj
    (
        "phi",
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phiIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating face flux field phi" << nl << endl;
    }
    else
    {
        Info<< "Reading face flux field phi" << nl << endl;
    }
    phi_.set
    (
        new surfaceScalarField
        (
            phiIOObj,
            linearInterpolate(rho_())*phiv_()
        )
    );

    if (!isStatic())
    {
        Info<< "Reading/calculating field Uf\n" << endl;
        Uf_.set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    runTime.timeName(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U_())
            )
        );
        fvc::makeRelative(phiv_(), U_());
    }

    pressureControl_.set
    (
        new pressureControl
        (
            p,
            rho_,
            solnControlPtr_->dict(),
            thermoPtr_->pRefValue()
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
            *thermoPtr_
        ).ptr()
    );
    turbulence_->store();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseSegregatedFlowSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {p_->name(), U_->name()};

    // Dependencies
    optionalDependencies.insert(p_->name(), {"fvMesh"});
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
}


bool Foam::fv::multiphaseSegregatedFlowSolver::initialise()
{
    if (!exists(phi_().objectPath()))
    {
        tmp<surfaceScalarField> rhof(fvc::interpolate(rho_()));
        this->fvOptions().makeRelative(rhof, phi_());
    }

    // Calculate initial Courant number
    maxCo_ =
        CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());

    return true;
}


scalar Foam::fv::multiphaseSegregatedFlowSolver::getMaxTimeStep()
{
    scalar maxMaxCo = 0;

    // Try region-specific maxCo, otherwise revert to global
    const dictionary& solnDict
    (
        solnControlPtr_->dict().found("maxCo")
      ? solnControlPtr_->dict()
      : mesh_.time().controlDict()
    );

    maxCoDataPtr_ = Function1<scalar>::New("maxCo", solnDict);

    maxMaxCo = maxCoDataPtr_->value(mesh_.time().timeIndex());

    return maxMaxCo/stabilise(maxCo_, SMALL)*mesh_.time().deltaTValue();
}


bool Foam::fv::multiphaseSegregatedFlowSolver::isFinalCorrector
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
        label nCorr = refCast<pimpleControl>(*solnControlPtr_).nCorrPISO();
        return (corrector >= nCorr-1);
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return false;
    }
}


const Foam::volScalarField&
Foam::fv::multiphaseSegregatedFlowSolver::buoyantRho() const
{
    return rho_();
}


void Foam::fv::multiphaseSegregatedFlowSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        implicitUCorrector_ =
            solnControlPtr_->dict().lookupOrDefault
            (
                "implicitUCorrector", false
            );
        implicitUSolve_ =
            solnControlPtr_->dict().lookupOrDefault
            (
                "implicitUSolve", implicitUCorrector_
            );
    }
    else if (correctorName == solverObject::outerCorrectorName)
    {
        if (!isStatic())
        {
            // Store divrhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
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


void Foam::fv::multiphaseSegregatedFlowSolver::assembleUEqnLHS()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();
    rho = thermoPtr_->rho();

    fv::options& fvOptions = this->fvOptions();

    adjustSplitBCs(phi_(), U_(), U);

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

    UEqnLHS_.ref() += fvm::ddt(rho, U);

    UEqnLHS_->relax();

    fvOptions.constrain(UEqnLHS_.ref());
}


void Foam::fv::multiphaseSegregatedFlowSolver::assemblepEqn()
{
    volVectorField& U = U_();
    surfaceScalarField& phiv = phiv_();
    volScalarField& rho = rho_();
    volScalarField& p = *p_;
    const volScalarField& psi = thermoPtr_->psi();
    fv::options& fvOptions = this->fvOptions();

    if (corr_ == 0 && nonOrthCorr_ == 0)
    {
        assembleUEqnLHS();
    }

    if (nonOrthCorr_ == 0)
    {
        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution
        p0_ = tmp<volScalarField>(new volScalarField(p));

        fvVectorMatrix& UEqnLHS(UEqnLHS_.ref());

        {
            tmp<volScalarField> A(UEqnLHS.A());
            tmp<volScalarField> AmH1(-UEqnLHS_().H1()+A());
            // Stabilise the case where H1 ~= -A
            rAU_ = new volScalarField("(1|At(U))", 1/max(AmH1, 0.01*A));
        }

        adjustSplitBCs(phi_(), U_(), p);

        const volScalarField& rAU = rAU_();
        HbyA_ =
            new volVectorField
            (
                constrainHbyA(rAU*(UEqnLHS.H() - UEqnLHS.H1()*U), U, p)
            );
        const volVectorField& HbyA = HbyA_();

        phivHbyA_ = new surfaceScalarField("phivHbyA", fvc::flux(HbyA));
        surfaceScalarField& phivHbyA(phivHbyA_.ref());
        rAUf_ = new surfaceScalarField("rAUf", fvc::interpolate(rAU));

        if (buoyant_)
        {
            phivg_ =
                rAUf_()
               *(
                    this->faceBuoyancyForce()
                  + thermoPtr_->surfaceTensionForce(U)
                )*mesh_.magSf();
        }
        else
        {
            phivg_ = rAUf_()*thermoPtr_->surfaceTensionForce(U)*mesh_.magSf();
        }

        if (Uf_.valid())
        {
            phivHbyA += fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, Uf_()); // Dynamic mesh
        }
        else
        {
            phivHbyA += fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phiv); // Static mesh
        }
        fvOptions.makeRelative(phivHbyA);

        closedVolume_ = p.needReference();

        phivHbyA += phivg_();

        // Update the fixedFluxPressure BCs to ensure flux consistency
        constrainPressure(p, U, phivHbyA, rAUf_(), fvOptions);

        pEqnComps_.setSize(thermoPtr_->thermos().size());

        forAll(thermoPtr_->thermos(), phasei)
        {
            const rhoThermo& thermoi =
                refCast<const rhoThermo>(thermoPtr_->thermos()[phasei]);
            tmp<volScalarField> trhoi = thermoi.rho();
            const volScalarField& rhoi = trhoi();

            if (!thermoi.isochoric())
            {
                pEqnComps_.set
                (
                    phasei,
                    new fvScalarMatrix
                    (
                        p, rho.dimensions()/dimTime*dimVolume
                    )
                );
                pEqnComps_[phasei] +=
                    fvc::ddt(rhoi)
                  + fvc::div(phiv, rhoi) - fvc::Sp(fvc::div(phiv), rhoi);
                // Compensate the div(meshPhi) in ddt(rhoi), since it cancels
                // in the two div terms above and the incompressible part of the
                // equation uses absolute fluxes
                if (!isStatic())
                {
                    pEqnComps_[phasei] -=
                        fvc::div
                        (
                            mesh_.phi(),
                            rhoi,
                            "div("+phiv.name()+','+rhoi.name()+')'
                        );
                }

                if (!thermoi.incompressible())
                {
                    pEqnComps_[phasei] += thermoi.psi()*correction(fvm::ddt(p));
                }
            }
        }

        if (!thermoPtr_->incompressible())
        {
            compressibility_() = fvc::domainIntegrate(psi);
        }

        pDDtEqn_ =
            new fvScalarMatrix
            (
              - (1/rho)*fvOptions(psi, p, rho.name())
              + fvc::div(phivHbyA)
            );
    }

    // Solve pressure
    pEqnIncomp_ =
        new fvScalarMatrix
        (
            pDDtEqn_()
          - fvm::laplacian(rAUf_, p)
        );

    pEqn_ = new fvScalarMatrix(pEqnIncomp_());

    forAll(thermoPtr_->thermos(), phasei)
    {
        if (pEqnComps_.set(phasei))
        {
            tmp<fvScalarMatrix> hmm
            (
                (
                    max(thermoPtr_->alphas()[phasei], scalar(0))
                   /thermoPtr_->thermos()[phasei].rho()
                )
               *pEqnComps_[phasei]
            );

            pEqn_.ref() += hmm;
        }
    }

    // Do not set reference in (transient), compressible case as there is a
    // dp/dt term
    if (p.needReference() && thermoPtr_->incompressible())
    {
        pEqn_->setReference
        (
            pRefCell_,
            pRefValue_
        );
    }
}


tmp<fvScalarMatrix>
Foam::fv::multiphaseSegregatedFlowSolver::assembleScalarMatrix
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
                finalIter_[solverObject::outerCorrectorName]
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


tmp<fvVectorMatrix>
Foam::fv::multiphaseSegregatedFlowSolver::assembleVectorMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if (fieldName == U_->name())
    {
        tmp<surfaceScalarField> tsnGrad;
        if (buoyant_)
        {
            tsnGrad = -this->faceBuoyancyForce(true);
        }
        else
        {
            tsnGrad = fvc::snGrad(*p_);
        }
        tsnGrad.ref() -= thermoPtr_->surfaceTensionForce(U_());
        if
        (
            (implicitUSolve_ && corr_ == 0)
         || (implicitUCorrector_ && corr_ != 0)
        )
        {
            return UEqnLHS_() + fvc::reconstruct(tsnGrad*mesh_.magSf());
        }
        else
        {
            tmp<fvVectorMatrix> UEqn
            (
                UEqnLHS_() + fvc::reconstruct(tsnGrad*mesh_.magSf())
            );
            U_() = UEqn().H()/UEqn().A();
            U_->correctBoundaryConditions();
            return tmp<fvVectorMatrix>();
        }
    }
    else
    {
        return tmp<fvVectorMatrix>();
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();
    if (solveName == p_->name())
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

            forAll(thermoPtr_->thermos(), phasei)
            {
                if (pEqnComps_.set(phasei))
                {
                    thermoPtr_->dgdts()[phasei] =
                        pos0(thermoPtr_->alphas()[phasei])
                       *(pEqnComps_[phasei] & p)
                       /thermoPtr_->thermos()[phasei].rho();
                }
                else
                {
                    thermoPtr_->dgdts()[phasei] *= 0;
                }
            }

            // Calculate the conservative flux
            phiv_() = phivHbyA_ + pEqnIncomp_->flux();
            if (!isStatic())
            {
                fvc::makeRelative(phiv_(), *U_);
            }

            pressureControl_->limit(p);

            if (closedVolume_)
            {
                if (thermoPtr_->incompressible())
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
                                thermoPtr_->rho() + thermoPtr_->psi()*(p-p0_())
                            )
                        )/compressibility_();
                }
            }

            if (!thermoPtr_->incompressible())
            {
                thermoPtr_->correctRho(p - p0_()); // Todo: make single phase thermo also work with dp instead of drho to be consistent with this
            }

            rho = thermoPtr_->rho();
        }
    }
    else if (solveName == U_->name())
    {
        // U corrector
        volVectorField& U = U_();

        pEqn_.clear();

        fvOptions.correct(U_());

        if (verbose_)
        {
            Info<< "mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }

        if (!isStatic())
        {
            Uf_() = fvc::interpolate(U);
            surfaceVectorField n(mesh_.Sf()/mesh_.magSf());
            Uf_() +=
                n
               *(
                    fvc::absolute(phiv_(), U)/mesh_.magSf()
                  - (n & Uf_())
                );
        }
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::solveContinuity()
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


void Foam::fv::multiphaseSegregatedFlowSolver::printIncoContinuityErrors
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


void Foam::fv::multiphaseSegregatedFlowSolver::printCompContinuityErrors
(
    const word& regionName
)
{
    const volScalarField& rho = rho_();

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr, globalContErr;

    sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rho - thermoPtr_->rho()))/totalMass
    ).value();

    globalContErr =
    (
        fvc::domainIntegrate(rho - thermoPtr_->rho())/totalMass
    ).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::multiphaseSegregatedFlowSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName && finalIter)
    {
        Info<< endl;
        // Calculate Courant for next time step
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), geometricOneField(), phiv_());
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::correctPhi(bool defaultCorrectPhi)
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
        phiv_() = mesh_.Sf() & Uf_();

        CorrectPhi
        (
            U_(),
            phi_(),
            *p_,
            rho_(),
            thermoPtr_->psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            divRhoU_(),
            *solnControlPtr_
        );

        // Make the fluxes relative to the mesh-motion
        fvc::makeRelative(phiv_(), U_());
    }
}


bool Foam::fv::multiphaseSegregatedFlowSolver::movePoints()
{
    calculateghFields(buoyant_);
    return true;
}


void Foam::fv::multiphaseSegregatedFlowSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    correctPhi(true);
    calculateghFields(buoyant_);
}


// ************************************************************************* //
