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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2020-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "foamCoupledControl.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamCoupledControl, 0);

    addToRunTimeSelectionTable
    (
        solutionControl,
        foamCoupledControl,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::foamCoupledControl::readControls()
{
    solutionControl::readControls(false);

    const dictionary& solverDict = dict();

    solveFlow_ = solverDict.lookupOrDefault<Switch>("solveFlow", true);
    nCorr_ = solverDict.lookupOrDefault<label>("nOuterCorrectors", 1);
    turbOnFinalIterOnly_ =
        solverDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", false);
    stabilityMode_ = solverDict.lookupOrDefault<label>("stabilityMode", 1);
    choiCorrection_ =
        solverDict.lookupOrDefault<Switch>("choiCorrection", false);
    procRelax_ = solverDict.lookupOrDefault<scalar>("procRelax", 1.0);

    meshQualityRelax_ =
        solverDict.lookupOrDefault<scalar>("meshQualityRelax", 1.0);

    localCellProcRelax_ =
        solverDict.lookupOrDefault<scalar>("localCellProcRelax", 0.7);
    if (meshQualityRelax_<1.0)
    {
        if (solverDict.found("meshQualityControls"))
        {
            const dictionary& meshQdict =
                solverDict.subDict("meshQualityControls");

            nonOrthoThreshold_ =
                meshQdict.lookupOrDefault<scalar>("nonOrthoThreshold", 60.0);
            skewnessThreshold_ =
                meshQdict.lookupOrDefault<scalar>("skewnessThreshold", 0.9);
        }
    }
}


bool Foam::foamCoupledControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.blockSolverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi = applyToField(variableName);
        if (fieldi != -1)
        {

            scalar residual = 0;
            const scalar firstResidual =
                maxResidual(iter(), residual);

            checked = true;

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = firstResidual;
            }

            const bool absCheck = residual < residualControl_[fieldi].absTol;
            bool relCheck = false;

            scalar relative = 0.0;
            if (!storeIni)
            {
                const scalar iniRes =
                    residualControl_[fieldi].initialResidual
                  + ROOTVSMALL;

                relative = residual/iniRes;
                relCheck = relative < residualControl_[fieldi].relTol;
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << variableName
                    << " outer iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldi].initialResidual
                    << ", abs tol = " << residual
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldi].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::foamCoupledControl::setFirstIterFlag
(
    const bool check,
    const bool force
)
{
    DebugInformation
        << "corr:" << corr_
        << endl;

    solutionControl::setFirstIterFlag(check && corr_ <= 1, force);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamCoupledControl::foamCoupledControl
(
    fvMesh& mesh,
    const word& dictName
)
:
    foamCoupledControl(mesh, mesh.thisDb(), dictName)
{
}


Foam::foamCoupledControl::foamCoupledControl
(
    fvMesh& mesh,
    const objectRegistry& obr,
    const word& dictName
)
:
    solutionControl(mesh, obr, dictName),
    solveFlow_(true),
    nCorr_(0),
    corr_(0),
    turbOnFinalIterOnly_(true),
    stabilityMode_(1),
    choiCorrection_(false),
    procRelax_(1.0),
    meshQualityRelax_(1.0),
    nonOrthoThreshold_(60.0),
    skewnessThreshold_(0.9),
    localCellProcRelax_(0.7),
    converged_(false)
{
    readControls();

    if (nCorr_ > 1)
    {
        Info<< nl;
        if (residualControl_.empty())
        {
            Info<< algorithmName_ << ": no residual control data found. "
                << "Calculations will employ " << nCorr_
                << " corrector loops" << nl << endl;
        }
        else
        {
            Info<< algorithmName_ << ": max iterations = " << nCorr_
                << endl;
            forAll(residualControl_, i)
            {
                Info<< "    field " << residualControl_[i].name << token::TAB
                    << ": relTol " << residualControl_[i].relTol
                    << ", tolerance " << residualControl_[i].absTol
                    << nl;
            }
            Info<< endl;
        }
    }
    else
    {
        Info<< nl << algorithmName_ << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamCoupledControl::~foamCoupledControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::foamCoupledControl::loop()
{
    readControls();

    corr_++;

    if (debug)
    {
        Info<< mesh_.time().timeName() << " loop: corr = " << corr_ << endl;
    }

    setFirstIterFlag();

    if (corr_ == nCorr_ + 1)
    {
        if ((!residualControl_.empty()) && (nCorr_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                 << nCorr_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                 << " iterations" << endl;

            mesh_.data::remove("finalIteration");

            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;

            mesh_.data::add("finalIteration", true);

            converged_ = true;
        }
    }
    else
    {
        if (corr_ <= nCorr_)
        {
            if (nCorr_!=1)
            {
                Info<< "Outer loop "
                     << mesh_.time().timeName()
                     << " | " << corr_ << endl;
            }

            if ((corr_ == nCorr_) && (nCorr_>1))
            {
                mesh_.data::add("finalIteration", true);
            }


            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
