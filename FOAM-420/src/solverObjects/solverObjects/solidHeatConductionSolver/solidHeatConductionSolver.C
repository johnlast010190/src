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

#include "solidHeatConductionSolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "sources/derived/GRFSource/GRFSource.H"
#include "sources/derived/solidRotationSource/solidRotationSource.H"
#include "finiteVolume/ddtSchemes/steadyStateDdtScheme/steadyStateDdtScheme.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidHeatConductionSolver, 0);
}
}

makeFvSolverOption(solidHeatConductionSolver);


// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //

void Foam::fv::solidHeatConductionSolver::updateAniAlpha()
{
    solidThermo& thermo = *thermoPtr_;

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    const coordinateSystem& coordinates = coorFramePtr_->coorSys();
    volVectorField KbyCp(thermo.Kappa()/cp);
    aniAlpha_->primitiveFieldRef() =
        coordinates.transformPrincipal(mesh_.C(), KbyCp);
    forAll(aniAlpha_->boundaryField(), patchi)
    {
        fvPatchField<symmTensor>& pf =
            aniAlpha_->boundaryFieldRef()[patchi];
        if (!isA<emptyFvPatch>(pf.patch()))
        {
            pf.forceAssign
            (
                coordinates.transformPrincipal
                (
                    mesh_.C().boundaryField()[patchi],
                    KbyCp.boundaryField()[patchi]
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidHeatConductionSolver::solidHeatConductionSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    coorFramePtr_(nullptr)
{
    const Time& runTime = mesh().time();

    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr);

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<solidThermo>
        (
            multiphaseThermo::lookupOrCreate(obr, phaseName_)
        );

    mesh_.schemes().setFluxRequired(thermoPtr_->he().name());
    mesh_.schemes().setFluxRequired(thermoPtr_->T().name());

    if (!thermoPtr_->isotropic())
    {
        Info<< "Adding reference frame for anisotropic solid\n" << endl;

        //TODO: We should make the anIso transport class
        // deal with the transformations itself and return the conductivity
        // tensor instead of the principal components

        word frameName = word::null;
        if (thermoPtr_->isDict("kappaModelCoeffs"))
        {
            frameName =
                thermoPtr_->subDict("kappaModelCoeffs").lookup<word>
                (
                    "referenceFrame"
                );
        }
        coorFramePtr_ = &coordinateFrame::New(mesh_, frameName);

        tmp<volVectorField> tkappaByCp = thermoPtr_->Kappa()/thermoPtr_->Cp();

        aniAlpha_.set
        (
            new volSymmTensorField
            (
                IOobject
                (
                    IOobject::groupName("Anialpha", phaseName_),
                    runTime.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedSymmTensor
                (
                    "zero",
                    tkappaByCp().dimensions(),
                    symmTensor::zero
                )
            )
        );

        updateAniAlpha();
    }

    IOobject betavIO
    (
        IOobject::groupName("betavSolid", phaseName_),
        runTime.timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (betavIO.typeHeaderOk<volScalarField>(true))
    {
        betav_.set(new volScalarField(betavIO, mesh_));
    }
    else
    {
        betavIO.readOpt() = IOobject::NO_READ;
        betav_.set
        (
            new volScalarField
            (
                betavIO,
                mesh_,
                dimensionedScalar("1", dimless, scalar(1.0))
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::solidHeatConductionSolver::initialise()
{
    // Initial diffusion number
    if (!isA<simpleControl>(*solnControlPtr_))
    {
        calcDiffusionNumber();
    }

    gradT_.set(fvc::grad(thermoPtr_->T()).ptr());

    // Suggest solidRotationSource instead of GRF
    fv::optionList& fvOptionList =
        dynamic_cast<fv::optionList&>(this->fvOptions());
    forAll(fvOptionList, optioni)
    {
        const fv::option& option = fvOptionList[optioni];
        if (isA<GRFSource>(option) && !isA<solidRotationSource>(option))
        {
            IOWarningInFunction(this->fvOptions())
                << "For solid regions, is is advisable to use a "
                << solidRotationSource::typeName
                << " instead of a " << GRFSource::typeName
                << " to model rotation."
                << nl << endl;
        }
    }

    return true;
}


void Foam::fv::solidHeatConductionSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    word hName(thermoPtr_->he().name());
    solveNames.append(hName);

    optionalDependencies.insert(hName, {"fvMesh"});

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
    correctorMembers.insert("nonOrthogonalCorrector:he", solveNames);

    requiredDependencies.insert("materialProperties", solveNames);
}


void Foam::fv::solidHeatConductionSolver::calcDiffusionNumber()
{
    solidThermo& thermo = *thermoPtr_;

    tmp<volScalarField> magKappa;
    if (thermo.isotropic())
    {
        magKappa = thermo.kappa();
    }
    else
    {
        magKappa = mag(thermo.Kappa());
    }

    maxDi_ =
        diffusionNumber
        (
            mesh_,
            mesh_.time(),
            thermo.rho(),
            magKappa()/thermo.Cp()
        );
}


scalar Foam::fv::solidHeatConductionSolver::getMaxTimeStep()
{
    scalar maxMaxDi = 0;
    // Try region-specific maxDi, otherwise revert to global
    if (solnControlPtr_->dict().found("maxDi"))
    {
        solnControlPtr_->dict().lookup("maxDi") >> maxMaxDi;
    }
    else
    {
        maxMaxDi = mesh_.time().controlDict().lookupOrDefault("maxDi", 10);
    }
    return maxMaxDi/maxDi_*mesh_.time().deltaTValue();
}


bool Foam::fv::solidHeatConductionSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if (correctorName == "nonOrthogonalCorrector:he")
    {
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return true;
    }
}


void Foam::fv::solidHeatConductionSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        // Update aniAlphas here so that they are available for neighbouring
        // regions

        if (!thermoPtr_->isotropic())
        {
            updateAniAlpha();
        }
    }
}


tmp<fvScalarMatrix>
Foam::fv::solidHeatConductionSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& solveName
)
{
    solidThermo& thermo = *thermoPtr_;
    if (fieldName == thermo.he().name())
    {
        fv::options& fvOptions = this->fvOptions();

        const volScalarField& rho = thermo.rho();

        volScalarField& he = thermo.he();
        const volScalarField& betav = betav_;

        word scheme = "laplacian(kappa";
        if (thermo.alpha().group() != word::null)
        {
            scheme += "." + thermo.alpha().group();
        }
        scheme += "," + thermo.T().name() + ")";

        tmp<fvScalarMatrix> thEqn;
        if (thermo.isotropic())
        {
            thEqn =
               -changeVariable
                (
                    fvm::laplacian(betav*thermo.kappa(), thermo.T(), scheme),
                    1/thermo.Cpv(),
                    he
                );
        }
        else
        {
            tmp<volScalarField> tCp = thermo.Cp();
            tmp<volSymmTensorField> betaKappa = betav*aniAlpha_()*tCp();
            tmp<volScalarField> rCpv = 1/tCp;
            if (he.member() == "e")
            {
                rCpv.ref() *= thermo.CpByCpv();
            }
            thEqn =
               -changeVariable
                (
                    fvm::laplacian(betaKappa, thermo.T(), scheme),
                    rCpv,
                    he
                );
        }

        const word ddtSchemeName = "ddt("+rho.name()+","+he.name()+')';
        if
        (
            !isA<steadyStateDdtScheme<scalar>>
            (
                ddtScheme<scalar>::New
                (
                    mesh_,
                    mesh_.schemes().ddtScheme(ddtSchemeName)
                )()
            )
        )
        {
            // Account for variation of rho in time derivative for non-
            // constant properties
            tmp<volScalarField> betavRho = betav*rho;
            volScalarField* betavRhoPtr = &(betavRho.ref());
            const volScalarField* rhoPtr = &rho;
            for (label nOldTime = 0; nOldTime < rho.nOldTimes(); nOldTime++)
            {
                betavRhoPtr = &betavRhoPtr->oldTime();
                rhoPtr = &rhoPtr->oldTime();
                *betavRhoPtr = betav*(*rhoPtr);
            }
            thEqn.ref() += fvm::ddt(betavRho, he, ddtSchemeName);
        }

        thEqn->relax();

        // Relaxation applied after GRF convection is added can slow down
        // conduction in the transverse direction, so we relax before adding
        // sources
        thEqn.ref() -= fvOptions(rho, he);

        fvOptions.constrain(thEqn.ref());

        return thEqn;
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }

}


void Foam::fv::solidHeatConductionSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (solveName == thermoPtr_->he().name())
    {
        this->fvOptions().correct(thermoPtr_->he());

        // Store the grad for use in non-orthog correction and BCs
        // Prevents duplicate calculation for anisotropic regions; we call it
        // regardless to allow consolidated group to contain both isotropic and
        // anisotropic regions
        gradT_.clear();
        gradT_.set(fvc::grad(thermoPtr_->T()).ptr());
    }
}


void Foam::fv::solidHeatConductionSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if
    (
        !isA<simpleControl>(*solnControlPtr_)
     && correctorName == solverObject::outerCorrectorName
     && finalIter
    )
    {
        // Calculate diffusion number for next time step
        calcDiffusionNumber();
    }
}

// ************************************************************************* //

