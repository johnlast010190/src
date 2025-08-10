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
    (c) 2019-2023 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "speciesConcentrationSolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "solidReactionThermo/solidReactionThermo.H"
#include "rhoReactionThermo/rhoReactionThermo.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(speciesConcentrationSolver, 0);
}
}

makeFvSolverOption(speciesConcentrationSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::speciesConcentrationSolver::speciesConcentrationSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    turbulencePtr_(nullptr),
    compositionPtr_(nullptr),
    inertIndex_(-1)
{
    Info<< "Reading thermophysical properties\n" << endl;
    //- Phasic and top-level thermo - will point to the same thing if phaseName_
    //  is null
    thermoPtr_ =
        &refCast<rhoThermo>
        (
            multiphaseThermo::lookupOrCreate(obr_, phaseName_)
        );
    globalThermoPtr_ =
        &refCast<rhoThermo>
        (
            basicThermo::lookupOrCreate(obr_, phaseName_)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::speciesConcentrationSolver::read(const dictionary& dict)
{
    // Read the Diffusion coefficients from thermo dict

    if (thermoPtr_ && compositionPtr_)
    {
        Di_.clear();
        Di_.resize(compositionPtr_->Y().size());
        if (thermoPtr_->phaseDict().isDict("speciesDiffusivity"))
        {
            const dictionary& diffDict =
                thermoPtr_->phaseDict().subDict("speciesDiffusivity");
            forAll(Di_, i)
            {
                if (diffDict.found(compositionPtr_->Y()[i].name()))
                {
                    Di_.set
                    (
                        i,
                        Function1<scalar>::New
                        (
                            compositionPtr_->Y()[i].name(),
                            diffDict
                        )
                    );
                }
            }
        }
    }
}


bool Foam::fv::speciesConcentrationSolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    if
    (
        obr_.foundObject<compressible::turbulenceModel>
        (
            compressible::turbulenceModel::propertiesName
        )
    )
    {
        turbulencePtr_ = &obr_.lookupObject<compressible::turbulenceModel>
        (
            compressible::turbulenceModel::propertiesName
        );
    }

    compositionPtr_ = &refCast<rhoReactionThermo>(*thermoPtr_).composition();

    word inertSpecies;
    if (thermoPtr_->phaseDict().found("inertSpecie"))
    {
        const word inertType
        (
            thermoPtr_->phaseDict().lookup<word>("inertSpecie")
        );
        inertSpecies = (inertType == "none") ? word::null : inertType;
    }
    else
    {
        const PtrList<volScalarField>& Y = compositionPtr_->Y();

        // Since the estimate of the inertSpecie doesn't have to be
        // exact the simple sum over field should be enough to compare
        // the fields.
        scalarField volSumY(Y.size());
        forAll(Y, i)
        {
            volSumY[i] =
                gSum(Y[i].primitiveField()) + gSum(Y[i].boundaryField());
        }
        inertSpecies = Y[findMax(volSumY)].member();
    }

    if (inertSpecies != word::null)
    {
        if (!compositionPtr_->species().found(inertSpecies))
        {
            FatalIOErrorInFunction(*thermoPtr_)
                << "Inert species " << inertSpecies
                << " not found in available species "
                << compositionPtr_->species()
                << exit(FatalIOError);
        }
        else
        {
            inertIndex_ = compositionPtr_->species()[inertSpecies];
        }
    }

    for (auto& Yi : compositionPtr_->Y())
    {
        multivariateConvectionFields_.add(Yi);

        // For the energy flux added in addSu
        mesh_.schemes().setFluxRequired(Yi.name());
    }

    read(dictionary(dict()));

    return true;
}


void Foam::fv::speciesConcentrationSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    const PtrList<volScalarField>& Y = compositionPtr_->Y();
    for (auto& Yi : Y)
    {
        solveNames.append(Yi.name());
        wordList deps;
        deps.append("fvMesh");
        word phiName = "phi";
        if (obr_.foundObject<surfaceScalarField>(phiName))
        {
            deps.append("U");
        }
        optionalDependencies.insert(Yi.name(), deps);

        correctorMembers.insert
        (
            "nonOrthogonalCorrector:"+Yi.member(), {Yi.name()}
        );
    }

    // Make the inert species dependent on all the others as we need them
    // to solve first
    if (inertIndex_ != -1)
    {
        DynamicList<word> activeSpecies;
        forAll(Y, speciesI)
        {
            if (speciesI != inertIndex_)
            {
                activeSpecies.append(Y[speciesI].name());
            }
        }
        requiredDependencies.set(Y[inertIndex_].name(), activeSpecies);
    }

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);

    requiredDependencies.insert("materialProperties", solveNames);
}


bool Foam::fv::speciesConcentrationSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if
    (
        correctorName.substr(0, correctorName.find(':'))
     == "nonOrthogonalCorrector"
    )
    {
        // Non-orthogonal correctors
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return false;
    }
}


tmp<fvScalarMatrix>
Foam::fv::speciesConcentrationSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    volScalarField& rho =
        obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));

    // If phi not found (no flow solver running), we include no convection
    tmp<fv::convectionScheme<scalar>> mvConvection;
    const surfaceScalarField* phiPtr = nullptr;
    word phiName = "phi";
    word YiName = addPhaseName("Yi");
    if (obr_.foundObject<surfaceScalarField>(phiName))
    {
        phiPtr = &obr_.lookupObjectRef<surfaceScalarField>(phiName);

        mvConvection =
            tmp<fv::convectionScheme<scalar>>
            (
                fv::convectionScheme<scalar>::New
                (
                    mesh_,
                    multivariateConvectionFields_,
                    *phiPtr,
                    mesh_.schemes().divScheme("div("+phiName+","+YiName+")")
                )
            );
    }

    PtrList<volScalarField>& Y = compositionPtr_->Y();
    forAll(Y, specieI)
    {
        if (Y[specieI].name() == fieldName)
        {
            if (compositionPtr_->active(specieI))
            {
                volScalarField& Yi = Y[specieI];

                fv::options& fvOptions = this->fvOptions();

                // Store diffusion coefficient for use by coupling BCs
                const word rhoDName = "rhoD_"+Yi.name();
                if (!obr_.foundObject<volScalarField>(rhoDName))
                {
                    tmp<volScalarField> rhoD
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                rhoDName,
                                mesh_.time().timeName(),
                                obr_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh_,
                            dimDensity*sqr(dimLength)/dimTime
                        )
                    );
                    regIOobject::store(rhoD.ptr());
                }
                volScalarField& rhoD =
                    obr_.lookupObjectRef<volScalarField>(rhoDName);

                // Laminar diffusion coefficient
                // Not (yet) integrated into thermo library as this depends on
                // the diffusion model (Fick vs multi-component, etc)
                if (Di_.set(specieI))
                {
                    // Specified diffusion vs temp
                    volScalarField D
                    (
                        IOobject
                        (
                            "D_"+YiName,
                            mesh_.time().timeName(),
                            obr_
                        ),
                        mesh_,
                        sqr(dimLength)/dimTime,
                        zeroGradientFvPatchScalarField::typeName
                    );
                    D.primitiveFieldRef() =
                        Di_[specieI].value(thermoPtr_->T().primitiveField())();
                    D.correctBoundaryConditions();

                    rhoD = rho*D;
                }
                else
                {
                    // If not specified, default to laminar viscosity (Reynolds
                    // analogy)
                    rhoD = thermoPtr_->mu();
                }

                // Turbulent contribution
                if (turbulencePtr_)
                {
                    rhoD += turbulencePtr_->Dmt();
                }

                // Tortuosity adjustment
                if (thermoPtr_->found("tortuosity"))
                {
                    scalar tortuosity =
                        readScalar(thermoPtr_->lookup("tortuosity"));
                    rhoD /= tortuosity;
                }

                if (specieI != inertIndex_)
                {
                    tmp<fvScalarMatrix> tYiEqn
                    (
                        new fvScalarMatrix
                        (
                            fvm::ddt(rho, Yi, "ddt("+rho.name()+",Yi)")
                          - fvm::laplacian
                            (
                                rhoD, Yi, "laplacian(rhoD_Yi,Yi)"
                            )
                        ==
                            fvOptions(rho, Yi)
                        )
                    );
                    if (mvConvection.valid())
                    {
                        tYiEqn.ref() += mvConvection->fvmDiv(*phiPtr, Yi);
                    }

                    if (thermoPtr_->found("porosity"))
                    {
                        scalar eps = readScalar(thermoPtr_->lookup("porosity"));
                        tYiEqn.ref() *= eps;
                    }

                    tYiEqn.ref().relax(YiName);
                    fvOptions.constrain(tYiEqn.ref());

                    dictName = YiName;
                    return tYiEqn;
                }
            }
        }
    }
    return tmp<fvScalarMatrix>();
}


void Foam::fv::speciesConcentrationSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    forAll(Y, specieI)
    {
        if (Y[specieI].name() == solveName)
        {
            if (specieI != inertIndex_ && compositionPtr_->active(specieI))
            {
                this->fvOptions().correct(Y[specieI]);
                Y[specieI].max(0.0);
            }
            else if (specieI == inertIndex_)
            {
                Y[inertIndex_] = scalar(1);
                forAll(Y, i)
                {
                    if (i != inertIndex_)
                    {
                        Y[inertIndex_] -= Y[i];
                    }
                }
                Y[inertIndex_].max(0.0);

                // Inert species is the last to be solved so we wait until here
                // to apply the density correction
                globalThermoPtr_->correct();

                // For stability, keep the solver density updated
                volScalarField& rho =
                    obr_.lookupObjectRef<volScalarField>("rho");
                rho = globalThermoPtr_->rho();
            }
            break;
        }
    }
}


void Foam::fv::speciesConcentrationSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{}


void Foam::fv::speciesConcentrationSolver::getSourceGraph
(
    wordList& fieldNames,
    HashTable<wordList>& sourceDependencies
)
{
    // We act as a source for the energy equation, to add the energy flux
    // due to mass transfer

    // Equations for which we are providing a source term
    fieldNames = {globalThermoPtr_->he().name()};

    // Require the mass fractions to have been solved before these terms are
    // added to energy equation
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    DynamicList<word> speciesNames;
    forAll(Y, specieI)
    {
        speciesNames.append(Y[specieI].name());
    }
    sourceDependencies.insert(fieldNames[0], speciesNames);
}


void Foam::fv::speciesConcentrationSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Adding terms to the energy equation to account for energy flux
    // due to mass transfer
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    const volScalarField& p = globalThermoPtr_->p();
    const volScalarField& T = globalThermoPtr_->T();
    const volScalarField& he = globalThermoPtr_->he();
    volScalarField hej("hej", he);

    forAll(Y, specieI)
    {
        if (compositionPtr_->active(specieI))
        {
            volScalarField& Yi = Y[specieI];

            const word rhoDName = "rhoD_"+Yi.name();
            volScalarField& rhoD =
                obr_.lookupObjectRef<volScalarField>(rhoDName);

            if (obr_.foundObject<objectRegistry>("materialModels"))
            {
                const objectRegistry& matObr = obr_.subRegistry("materialModels");
                const materialTables& matTable =
                    matObr.lookupObject<materialTables>("materialTables");
                hej.forceAssign(matTable(HEModel::typeName, word::null, Yi.name())());
            }
            else
            {
                // Get the enthalpy or energy for this species
                forAll(thermoPtr_->T(), celli)
                {
                    hej[celli] =
                        compositionPtr_->HE
                        (
                            specieI, p[celli]+thermoPtr_->pRefValue(), T[celli]
                        );
                }
                volScalarField::Boundary& hejBf = hej.boundaryFieldRef();
                forAll(hejBf, patchi)
                {
                    scalarField& hejp = hejBf[patchi];
                    const scalarField& pp = p.boundaryField()[patchi];
                    const scalarField& Tp = T.boundaryField()[patchi];

                    forAll(hejp, facei)
                    {
                        hejp[facei] =
                            compositionPtr_->HE
                            (
                                specieI,
                                pp[facei]+thermoPtr_->pRefValue(),
                                Tp[facei]
                            );
                    }
                }
            }

            // Species diffusive mass flux
            tmp<surfaceScalarField> tJi =
               -fvm::laplacian
                (
                    rhoD, Yi, "laplacian(rhoD_Yi,Yi)"
                )->flux();
            const surfaceScalarField& Ji = tJi();

            tmp<fv::convectionScheme<scalar>> convScheme =
                fv::convectionScheme<scalar>::New
                (
                    mesh_,
                    Ji,
                    mesh_.schemes().divScheme
                    (
                        "div(Ji,"
                        +IOobject::groupName
                        (
                            eqn.psi().member()+"i",
                            eqn.psi().group()
                        )
                        +")"
                    )
                );

            // We want convection of the species energy, div(Ji, hej), but
            // treating implicitly through he.
            // eqn is RHS, so subtract
            eqn -=
                convScheme->fvmDiv
                (
                    Ji
                   *convScheme->interpolate(Ji, hej)
                   /convScheme->interpolate(Ji, eqn.psi()),
                    eqn.psi()
                );
        }
    }
}


// ************************************************************************* //
