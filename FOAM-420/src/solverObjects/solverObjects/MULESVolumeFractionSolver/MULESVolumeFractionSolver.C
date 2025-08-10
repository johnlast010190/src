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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2022-2024 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "MULESVolumeFractionSolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "solidReactionThermo/solidReactionThermo.H"
#include "rhoReactionThermo/rhoReactionThermo.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "fvMatrices/solvers/MULES/CMULES.H"
#include "algorithms/subCycle/subCycle.H"
#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(MULESVolumeFractionSolver, 0);
}
}

makeFvSolverOption(MULESVolumeFractionSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MULESVolumeFractionSolver::MULESVolumeFractionSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    thermoPtr_(nullptr),
    passiveIndex_(-1),
    phivPtr_(nullptr),
    phiPtr_(nullptr),
    maxAlphaCo_(0),
    solnControlPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::MULESVolumeFractionSolver::read(const dictionary& dict)
{}


bool Foam::fv::MULESVolumeFractionSolver::initialise()
{
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr_));

    passiveIndex_ = thermoPtr_->fractions().passiveIndex();

    read(dictionary(dict()));

    phivPtr_ = &mesh_.lookupObject<surfaceScalarField>("phiv");
    phiPtr_ = &mesh_.lookupObjectRef<surfaceScalarField>("phi");
    UPtr_ = &mesh_.lookupObject<volVectorField>("U");

    updateTimeSchemeInfo();

    return true;
}


void Foam::fv::MULESVolumeFractionSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    const wordList& phases = thermoPtr_->phases();
    solveNames.append("alphas");
    derivedFields.insert("alphas", phases);
    optionalDependencies.insert("alphas", {"fvMesh"});
    optionalDependencies.insert("p", {"alphas"});
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


scalar Foam::fv::MULESVolumeFractionSolver::getMaxTimeStep()
{
    if (!LTS_)
    {
        // Try region-specific maxCo, otherwise revert to global
        const dictionary& maxCoDict
        (
            solnControlPtr_->dict().found("maxAlphaCo")
          ? solnControlPtr_->dict()
          : mesh_.time().controlDict()
        );
        maxAlphaCoDataPtr_ = Function1<scalar>::New("maxAlphaCo", maxCoDict);

        scalar maxMaxAlphaCo =
            maxAlphaCoDataPtr_->value(mesh_.time().timeIndex());
        return
            maxMaxAlphaCo/stabilise(maxAlphaCo_, SMALL)
           *mesh_.time().deltaTValue();
    }
    else
    {
        return GREAT;
    }
}


bool Foam::fv::MULESVolumeFractionSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    return true;
}


void Foam::fv::MULESVolumeFractionSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        const dictionary& alphaControls = mesh().solution().solverDict("alpha");
        nAlphaSubCycles_ = readLabel(alphaControls.lookup("nAlphaSubCycles"));
        cAlpha_ = readScalar(alphaControls.lookup("cAlpha"));
        nAlphaCorr_ = alphaControls.lookupOrDefault<label>("nAlphaCorr", 1);
        MULESCorr_ = alphaControls.lookupOrDefault<Switch>("MULESCorr", false);
        alphaApplyPrevCorr_ =
            alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false);
        if (!alphaApplyPrevCorr_)
        {
            alphaPhiCorr0_.clear();
        }

        updateTimeSchemeInfo();

        if (MULESCorr_)
        {
            forAll(thermoPtr_->alphas(), phasei)
            {
                if (phasei != passiveIndex_ || thermoPtr_->alphas().size() > 2)
                {
                    mesh().schemes().setFluxRequired
                    (
                        thermoPtr_->alphas()[phasei].name()
                    );
                }
            }
        }

    }
}


void Foam::fv::MULESVolumeFractionSolver::updateTimeSchemeInfo()
{
    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh_,
                mesh_.schemes().ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            crankNicolson_ = false;
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            crankNicolson_ = true;
            if (mesh_.time().subCycling())
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            ocCoeff =
                refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                .ocCoeff();
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    if (!crankNicolson_)
    {
        alphaPhis_.clear();
    }
    else
    {
        if (!alphaPhis_.size())
        {
            alphaPhis_.setSize(thermoPtr_->phases().size());
            forAll(alphaPhis_, phasei)
            {
                alphaPhis_.set
                (
                    phasei,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "alphaPhi."+thermoPtr_->phases()[phasei],
                            mesh_.time().timeName(),
                            obr_,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        (*phivPtr_)
                       *fvc::interpolate(thermoPtr_->alphas()[phasei])
                    )
                );
            }
        }
    }

    // Set the time blending factor, 1 for Euler
    cnCoeff_ = 1.0/(1.0 + ocCoeff);

    LTS_ = fv::localEulerDdt::enabled(mesh_);

    if (LTS_ && !rDeltaT_.valid())
    {
        Info<< "Using LTS" << endl;

        rDeltaT_.set
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("one", dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
    else if (!LTS_)
    {
        rDeltaT_.clear();
    }
}

void Foam::fv::MULESVolumeFractionSolver::calcSuSp
(
    const label phasei,
    autoPtr<volScalarField::Internal>& pSu,
    autoPtr<volScalarField::Internal>& pSp
)
{
    const volScalarField& alpha = thermoPtr_->alphas()[phasei];
    const volScalarField& dgdt = thermoPtr_->dgdts()[phasei];

    pSu.set
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimVelocity/dimLength
        )
    );

    pSp.set
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dgdt.dimensions()
        )
    );

    volScalarField::Internal& Su(pSu());
    volScalarField::Internal& Sp(pSp());

    forAll(dgdt, celli)
    {
        if (dgdt[celli] < 0.0 && alpha[celli] > 0.0)
        {
            Sp[celli] =  dgdt[celli]*alpha[celli];
            Su[celli] = -dgdt[celli]*alpha[celli];
        }
        else if (dgdt[celli] > 0.0 && alpha[celli] < 1.0)
        {
            Sp[celli] = -dgdt[celli]*(1.0 - alpha[celli]);
            Su[celli] = 0.0;
        }
        else
        {
            Sp[celli] = 0.0;
            Su[celli] = 0.0;
        }
    }

    forAll(thermoPtr_->phases(), phase2i)
    {
        const volScalarField& alpha2 = thermoPtr_->alphas()[phase2i];

        if (&alpha2 != &alpha)
        {
            const scalarField& dgdt2 = thermoPtr_->dgdts()[phase2i];

            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*alpha2[celli];
                }
            }
        }
    }
}


void Foam::fv::MULESVolumeFractionSolver::solveAlphas()
{
    word alphaScheme("div(phiv,alpha)");
    word alpharScheme("div(phivrb,alpha)");

    const surfaceScalarField& phiv(*phivPtr_);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mag(phiv/mesh_.magSf()));
    phic = min(cAlpha_*phic, max(phic));

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();
    forAll(phicBf, patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp.forceAssign(0);
        }
    }

    // Calculate the Crank-Nicolson off-centred volumetric flux
    tmp<surfaceScalarField> phivCN(phiv);
    if (crankNicolson_)
    {
        phivCN = cnCoeff_*phiv + (1.0 - cnCoeff_)*phiv.oldTime();
        phivCN->rename(phiv.name());
    }

    const wordList& phases = thermoPtr_->phases();

    volScalarField::Internal divU(fvc::div(phiv));

    if (MULESCorr_)
    {
        volScalarField sumAlpha
        (
            "sumAlpha",
            phiv.db(),
            mesh_,
            dimless,
            0
        );

        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_ && phases.size() <= 2)
            {
                continue;
            }

            volScalarField& alpha = thermoPtr_->alphas()[phasei];

            autoPtr<volScalarField::Internal> Su, Sp;
            calcSuSp(phasei, Su, Sp);

            fvScalarMatrix alphaEqn
            (
                (
                    LTS_
                  ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                  : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                )
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh_,
                    phivCN,
                    upwind<scalar>(mesh_, phivCN)
                ).fvmDiv(phivCN, alpha)
             ==
                Su() + fvm::Sp(Sp() + divU, alpha)
            );

            alphaEqn.solve();

            tmp<surfaceScalarField> talphaPhiUD(alphaEqn.flux());
            if (crankNicolson_)
            {
                alphaPhis_[phasei] = talphaPhiUD();
            }
            else
            {
                alphaPhis_.setSize(thermoPtr_->phases().size());
                alphaPhis_.set(phasei, new surfaceScalarField(talphaPhiUD()));
            }

            if
            (
                alphaApplyPrevCorr_
             && alphaPhiCorr0_.size()
             && alphaPhiCorr0_.set(phasei)
            )
            {
                MULES::correct
                (
                    alpha,
                    alphaPhis_[phasei],
                    alphaPhiCorr0_[phasei],
                    1,
                    0
                );

                alphaPhis_[phasei] += alphaPhiCorr0_[phasei];
            }

            if (alphaApplyPrevCorr_)
            {
                // Cache the upwind-flux
                if (!alphaPhiCorr0_.size())
                {
                    alphaPhiCorr0_.setSize(thermoPtr_->phases().size());
                }
                alphaPhiCorr0_.set(phasei, talphaPhiUD.ptr());
            }

            if (passiveIndex_ != -1)
            {
                sumAlpha += alpha;
            }
        }

        if (passiveIndex_ != -1)
        {
            volScalarField& alpha = thermoPtr_->alphas()[passiveIndex_];
            alpha = 1-sumAlpha;
        }
    }

    PtrList<surfaceScalarField> alphaPhiCorrs(phases.size());
    PtrList<surfaceScalarField> alphaPhiUns(phases.size());

    for (label aCorr = 0; aCorr < nAlphaCorr_; aCorr++)
    {
        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_ && phases.size() <= 2)
            {
                continue;
            }

            const volScalarField& alpha = thermoPtr_->alphas()[phasei];

            tmp<surfaceScalarField> alphaPhiUn
            (
                new surfaceScalarField
                (
                    IOobject::groupName(phiv.name(), alpha.group()),
                    fvc::flux
                    (
                        phivCN(),
                        crankNicolson_
                      ? cnCoeff_*alpha + (1.0 - cnCoeff_)*alpha.oldTime()
                      : tmp<volScalarField>(alpha),
                        alphaScheme
                    )
                )
            );

            forAll(phases, phase2i)
            {
                const volScalarField& alpha2 = thermoPtr_->alphas()[phase2i];

                if (&alpha2 != &alpha)
                {
                    surfaceScalarField phir
                    (
                        phic*thermoPtr_->nHatf(alpha, alpha2)
                    );

                    alphaPhiUn.ref() +=
                        fvc::flux
                        (
                            -fvc::flux(-phir, alpha2, alpharScheme),
                            alpha,
                            alpharScheme
                        );
                }
            }

            if (MULESCorr_)
            {
                alphaPhiCorrs.set
                (
                    phasei,
                    alphaPhiUn() - alphaPhis_[phasei]
                );
                alphaPhiUns.set(phasei, alphaPhiUn);
                surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

                autoPtr<volScalarField::Internal> Su, Sp;
                calcSuSp(phasei, Su, Sp);

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);
                    MULES::limitCorr
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhiUns[phasei],
                        alphaPhiCorr,
                        Sp(),
                        (-Sp()*alpha)(),
                        1,
                        0
                    );
                }
                else
                {
                    MULES::limitCorr
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhiUns[phasei],
                        alphaPhiCorr,
                        Sp(),
                        (-Sp()*alpha)(),
                        1,
                        0
                    );
                }
            }
            else
            {
                alphaPhiCorrs.set(phasei, alphaPhiUn.ptr());
                surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    // NOTE: This converts alphaPhiCorr from a flux to a correction flux
                    MULES::limit
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        phiv,
                        alphaPhiCorr,
                        zeroField(),
                        zeroField(),
                        1,
                        0,
                        true
                    );
                }
                else
                {
                    // NOTE: This converts alphaPhiCorr from a flux to a correction flux
                    MULES::limit
                    (
                        1.0/mesh_.time().deltaT().value(),
                        geometricOneField(),
                        alpha,
                        phiv,
                        alphaPhiCorr,
                        zeroField(),
                        zeroField(),
                        1,
                        0,
                        true
                    );
                }
            }
        }

        if (passiveIndex_ < 0 || phases.size() > 2)
        {
            // We still need to limit the sum even when there is a passive phase
            // if there are more than 2 phases, as there could be interfaces
            // where the passive phase is absent and can't absorb a discrepancy

            MULES::limitSum(alphaPhiCorrs);
        }

        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_)
            {
                continue;
            }

            volScalarField& alpha = thermoPtr_->alphas()[phasei];

            autoPtr<volScalarField::Internal> Su, Sp;
            calcSuSp(phasei, Su, Sp);

            if (MULESCorr_)
            {
                volScalarField alpha0("alpha0", alpha);

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    MULES::correct
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhiUns[phasei],
                        alphaPhiCorrs[phasei],
                        Sp(),
                        (-Sp()*alpha)()
                    );
                }
                else
                {
                    MULES::correct
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhiUns[phasei],
                        alphaPhiCorrs[phasei],
                        Sp(),
                        (-Sp()*alpha)()
                    );
                }

                // Under-relax the correction for all but the 1st corrector
                if (aCorr == 0)
                {
                    alphaPhis_[phasei] += alphaPhiCorrs[phasei];
                }
                else
                {
                    alpha = 0.5*alpha + 0.5*alpha0;
                    alpha.correctBoundaryConditions();
                    alphaPhis_[phasei] += 0.5*alphaPhiCorrs[phasei];
                }

            }
            else
            {
                tmp<surfaceScalarField> talphaPhi =
                    upwind<scalar>(mesh_, phiv).flux(alpha)
                  + alphaPhiCorrs[phasei];
                if (crankNicolson_)
                {
                    alphaPhis_[phasei] = talphaPhi;
                }
                else
                {
                    if (!alphaPhis_.size())
                    {
                        alphaPhis_.setSize(thermoPtr_->phases().size());
                    }
                    alphaPhis_.set(phasei, talphaPhi);
                }

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    MULES::explicitSolve
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhis_[phasei],
                        Sp(),
                        // Divergence term is handled explicitly to be
                        // consistent with the explicit transport solution
                        (Su() + divU*min(alpha.internalField(), scalar(1)))()
                    );
                }
                else
                {
                    MULES::explicitSolve
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhis_[phasei],
                        Sp(),
                        // Divergence term is handled explicitly to be
                        // consistent with the explicit transport solution
                        (Su() + divU*min(alpha.internalField(), scalar(1)))()
                    );
                }
            }
        }
    }

    surfaceScalarField newPhi
    (
        "newPhi",
        phiPtr_->db(),
        mesh_,
        phiPtr_->dimensions(),
        0
    );
    volScalarField sumAlpha
    (
        "sumAlpha",
        phiv.db(),
        mesh_,
        dimless,
        0
    );
    autoPtr<surfaceScalarField> sumAlphaPhi;
    if (passiveIndex_ != -1)
    {
        sumAlphaPhi.set
        (
            new surfaceScalarField
            (
                "sumAlphaPhi",
                phiv.db(),
                mesh_,
                phiv.dimensions(),
                0
            )
        );
    }

    forAll(phases, phasei)
    {
        if (phasei == passiveIndex_)
        {
            continue;
        }

        volScalarField& alpha = thermoPtr_->alphas()[phasei];

        surfaceScalarField& alphaPhi = alphaPhis_[phasei];

        if (crankNicolson_)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi =
                (alphaPhi - (1.0 - cnCoeff_)*alphaPhi.oldTime())/cnCoeff_;
        }

        newPhi +=
            fvc::interpolate(thermoPtr_->thermos()[phasei].rho())*alphaPhi;

        sumAlpha += alpha;
        if (passiveIndex_ != -1)
        {
            sumAlphaPhi() += alphaPhi;
        }

        if (alphaApplyPrevCorr_)
        {
            alphaPhiCorr0_.set
            (
                phasei,
                alphaPhiCorrs[phasei] - alphaPhiCorr0_[phasei]
            );
            alphaPhiCorr0_[phasei].rename("alphaPhiCorr0");
        }

        Info<< alpha.group() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;
    }

    if (passiveIndex_ != -1)
    {
        volScalarField& alpha = thermoPtr_->alphas()[passiveIndex_];
        alpha = 1-sumAlpha;

        // Calculate the end-of-time-step mass flux (phiv rather than phiCN)
        newPhi +=
            fvc::interpolate(thermoPtr_->thermos()[passiveIndex_].rho())
            *(phiv-sumAlphaPhi());
    }
    else
    {
        Info<< "Phase-sum volume fraction, min, max = "
            << sumAlpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(sumAlpha).value()
            << ' ' << max(sumAlpha).value()
            << endl;
    }
    *phiPtr_ = newPhi;

    if (!crankNicolson_)
    {
        alphaPhis_.clear();
    }

    thermoPtr_->fractions().recomputeCombinedAlphas();
}


void Foam::fv::MULESVolumeFractionSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    const Time& runTime = mesh().time();

    volScalarField& alpha = thermoPtr_->alphas().first();

    if (nAlphaSubCycles_ > 1)
    {
        surfaceScalarField& phi(*phiPtr_);
        surfaceScalarField phiSum(0.0*phi);
        dimensionedScalar totalDeltaT = runTime.deltaT();

        tmp<volScalarField> trSubDeltaT;
        if (LTS_)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles_);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles_);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();
            phiSum += (runTime.deltaT()/totalDeltaT)*phi;
        }

        phi = phiSum;
    }
    else
    {
        solveAlphas();
    }

    thermoPtr_->correct();
}


void Foam::fv::MULESVolumeFractionSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName && finalIter)
    {
        Info<< endl;

        if (!LTS_)
        {
            // Calculate Courant number at interface

            scalar maxAlphaCo_ = 0.0;
            scalar meanAlphaCoNum = 0.0;

            if (mesh().nInternalFaces())
            {
                scalarField sumPhi
                (
                    thermoPtr_->nearInterface()().primitiveField()
                   *fvc::surfaceSum(mag(*phivPtr_))().primitiveField()
                );

                maxAlphaCo_ =
                    0.5*gMax(sumPhi/mesh().V().field())*mesh().time().deltaTValue();

                meanAlphaCoNum =
                    0.5
                   *(
                        gSum(sumPhi)/gSum(mesh().V().field())
                    )*mesh().time().deltaTValue();
            }

            Info<< "Interface Courant Number mean: " << meanAlphaCoNum
                << " max: " << maxAlphaCo_ << endl;
        }
        else
        {
            const dictionary& pimpleDict = solnControlPtr_->dict();

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

            volScalarField rDeltaT0("rDeltaT0", rDeltaT_());

            // Set the reciprocal time-step from the local Courant number
            rDeltaT_->ref() =
                max
                (
                    1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
                    fvc::surfaceSum(mag(*phivPtr_))()()
                   /((2*maxCo)*mesh_.V())
                );

            // Create a compound 'interface' indicater that it the product
            // of alphas
            volScalarField alpha(thermoPtr_->alphas()[0]);
            for
            (
                label phasei = 1;
                phasei < thermoPtr_->alphas().size();
                phasei++
            )
            {
                alpha *= thermoPtr_->alphas()[phasei];
            }

            if (maxAlphaCo < maxCo)
            {
                // Further limit the reciprocal time-step
                // in the vicinity of the interface

                volScalarField alphaBar(fvc::average(alpha));

                rDeltaT_->ref() =
                    max
                    (
                        rDeltaT_(),
                        pos0(alphaBar() - alphaSpreadMin)
                       *pos0(alphaSpreadMax - alphaBar())
                       *fvc::surfaceSum(mag(*phivPtr_))()()
                       /((2*maxAlphaCo)*mesh_.V())
                    );
            }

            // Update the boundary values of the reciprocal time-step
            rDeltaT_().correctBoundaryConditions();

            Info<< "Flow time scale min/max = "
                << gMin(1/rDeltaT_().primitiveField())
                << ", " << gMax(1/rDeltaT_().primitiveField()) << endl;

            if (rDeltaTSmoothingCoeff < 1.0)
            {
                fvc::smooth(rDeltaT_(), rDeltaTSmoothingCoeff);
            }

            if (nAlphaSpreadIter > 0)
            {
                fvc::spread
                (
                    rDeltaT_(),
                    alpha,
                    nAlphaSpreadIter,
                    alphaSpreadDiff,
                    alphaSpreadMax,
                    alphaSpreadMin
                );
            }

            if (nAlphaSweepIter > 0)
            {
                fvc::sweep(rDeltaT_(), alpha, nAlphaSweepIter, alphaSpreadDiff);
            }

            Info<< "Smoothed flow time scale min/max = "
                << gMin(1/rDeltaT_().primitiveField())
                << ", " << gMax(1/rDeltaT_().primitiveField()) << endl;

            // Limit rate of change of time scale
            // - reduce as much as required
            // - only increase at a fraction of old time scale
            if
            (
                rDeltaTDampingCoeff < 1.0
             && mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1
            )
            {
                rDeltaT_() =
                    max
                    (
                        rDeltaT_(),
                        (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
                    );

                Info<< "Damped flow time scale min/max = "
                    << gMin(1/rDeltaT_().primitiveField())
                    << ", " << gMax(1/rDeltaT_().primitiveField()) << endl;
            }
        }
    }
}


void Foam::fv::MULESVolumeFractionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    alphaPhiCorr0_.clear();
}

// ************************************************************************* //
