/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeControl.H"
#include "db/dictionary/dictionary.H"
#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "solverObject/solverObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(runTimeControl, 0);
    addToRunTimeSelectionTable(functionObject, runTimeControl, dictionary);
}
}
}

Foam::Enum
<
    Foam::functionObjects::runTimeControls::runTimeControl::satisfiedAction
>
Foam::functionObjects::runTimeControls::runTimeControl::satisfiedActionNames
{
    { satisfiedAction::ABORT, "abort"},
    { satisfiedAction::END, "end"},
    { satisfiedAction::SET_TRIGGER, "setTrigger"},
    { satisfiedAction::MODIFY_SOLVER, "modifySolver"},
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeControl::runTimeControl
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    conditions_(),
    groupMap_(),
    nWriteStep_(0),
    writeStepI_(0),
    satisfiedAction_(satisfiedAction::END),
    triggerIndex_(labelMin),
    active_(getObjectProperty(name, "active", true)),
    functionsOnStopping_
    (
        dict.lookupOrDefault<Switch>("runFunctionObjects", false)
    ),
    solveInterval_(new Function1Types::Constant<scalar>("solveInterval", 1))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::runTimeControl::read
(
    const dictionary& dict
)
{
    if (functionObject::postProcess)
    {
        Info<< "Deactivated " << name()
            << " function object for post-processing"
            << endl;

        return false;
    }


    if (fvMeshFunctionObject::read(dict))
    {
        const dictionary& conditionsDict = dict.subDict("conditions");
        const wordList conditionNames(conditionsDict.toc());
        conditions_.setSize(conditionNames.size());

        label uniqueGroupi = 0;
        forAll(conditionNames, conditioni)
        {
            const word& conditionName = conditionNames[conditioni];
            const dictionary& dict = conditionsDict.subDict(conditionName);

            conditions_.set
            (
                conditioni,
                runTimeCondition::New(conditionName, obr_, dict, *this)
            );

            label groupi = conditions_[conditioni].groupID();

            if (groupMap_.insert(groupi, uniqueGroupi))
            {
                ++uniqueGroupi;
            }
        }

        dict.readIfPresent("nWriteStep", nWriteStep_);

        // Check that some conditions are set
        if (conditions_.empty())
        {
            Info<< type() << " " << name() << " output:" << nl
                << "    No conditions present" << nl
                << endl;
        }
        else
        {
            // Check that at least one condition is active
            bool check = false;
            for (const auto& condition : conditions_)
            {
                if (condition.active())
                {
                    check = true;
                    break;
                }
            }

            if (!check)
            {
                Info<< type() << " " << name() << " output:" << nl
                    << "    All conditions are inactive" << nl
                    << endl;
            }
        }

        // Set the action to perform when all conditions are satisfied
        // - set to end for backwards compatibility with v1806
        if (dict.found("satisfiedAction"))
        {
            satisfiedAction_ =
            satisfiedActionNames.read
            (
                dict.lookup("satisfiedAction")
            );
        }

        if (satisfiedAction_ == satisfiedAction::SET_TRIGGER)
        {
            triggerIndex_ = dict.lookup<label>("trigger");
        }
        else if (satisfiedAction_ == satisfiedAction::MODIFY_SOLVER)
        {
            if (
                !obr_.foundObject<solverObject>("flowSolver") and
                !obr_.foundObject<solverObject>("turbulenceSolver")
            )
            {
                satisfiedAction_ = satisfiedAction::END;
                Warning << "flowSolver and turbulenceSolver not found so satisfiedAction is forced to end" << endl;
            }
            else
            {
                const dictionary& modifySolverDict = dict.subDict("modifySolverCoeffs");
                solveInterval_ = Function1<scalar>::New("solveInterval", modifySolverDict);
            }
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::execute()
{
    if (!active_)
    {
        return true;
    }

    Info<< type() << " " << name() << " output:" << nl;

    // IDs of satisfied conditions
    DynamicList<label> IDs(conditions_.size());

    // Run stops only if all conditions within a group are satisfied
    List<bool> groupSatisfied(groupMap_.size(), true);
    List<bool> groupActive(groupMap_.size(), false);

    forAll(conditions_, conditioni)
    {
        runTimeCondition& condition = conditions_[conditioni];

        if (condition.active())
        {
            bool conditionSatisfied = condition.apply();

            label groupi = condition.groupID();

            auto conditionIter = groupMap_.cfind(groupi);

            if (!conditionIter.found())
            {
                FatalErrorInFunction
                    << "group " << groupi << " not found in map"
                    << Foam::abort(FatalError);
            }

            if (conditionSatisfied)
            {
                IDs.append(conditioni);

                groupActive[conditionIter()] = true;

                if (groupi == -1)
                {
                    // Condition not part of a group - only requires this to be
                    // satisfied for completion flag to be set
                    groupSatisfied[conditionIter()] = true;
                    break;
                }
            }
            else
            {
                groupSatisfied[conditionIter()] = false;
            }
        }
    }

    bool done = false;
    forAll(groupSatisfied, groupi)
    {
        if (groupSatisfied[groupi] && groupActive[groupi])
        {
            done = true;
            break;
        }
    }

    if (done)
    {
        for (label conditioni : IDs)
        {
            Info<< "    " << conditions_[conditioni].type() << ": "
                << conditions_[conditioni].name()
                << " condition satisfied" << nl;
        }

        switch (satisfiedAction_)
        {
            case satisfiedAction::ABORT:
            case satisfiedAction::END:
            {
                // Set to write a data dump or finalise the calculation
                Time& time = const_cast<Time&>(time_);

                if (writeStepI_ < nWriteStep_ - 1)
                {
                    ++writeStepI_;
                    Info<< "    Writing fields - step " << writeStepI_ << nl;
                    time.writeNow();
                }
                else
                {
                    Info<< "    Stopping calculation" << nl
                        << "    Writing fields";

                    if (nWriteStep_ != 0)
                    {
                        Info<< " - final step" << nl;
                    }
                    else
                    {
                        Info<< nl;
                    }

                    Info<< endl;
                    active_ = false;

                    // Write any registered objects and set the end-time
                    time.writeAndEnd();

                    if (functionsOnStopping_)
                    {
                        // Trigger any function objects
                        time.functionObjects().execute();
                    }

                    if (satisfiedAction_ == satisfiedAction::ABORT)
                    {
                        FatalErrorInFunction
                            << "Abort triggered"
                            << exit(FatalError);
                    }
                }
                break;
            }
            case satisfiedAction::SET_TRIGGER:
            {
                Info<< "    Setting trigger " << triggerIndex_ << nl;
                setTrigger(triggerIndex_);

                // Deactivate the model
                active_ = false;
                setProperty("active", active_);
                break;
            }
            case satisfiedAction::MODIFY_SOLVER:
            {
                Info<< "    Performing modifySolver to update solveInterval" << nl;

                if (solveInterval_->type() == "table")
                {
                    const Time& currentTime = obr_.time();
                    solveInterval_->shiftTimeBase(currentTime);
                }

                if (obr_.foundObject<solverObject>("flowSolver")) {
                    solverObject& flowSolver = obr_.lookupObjectRef<solverObject>("flowSolver");
                    flowSolver.resetSolveInterval(solveInterval_);
                }

                if (obr_.foundObject<solverObject>("turbulenceSolver")) {
                    solverObject& turbulenceSolver = obr_.lookupObjectRef<solverObject>("turbulenceSolver");
                    turbulenceSolver.resetSolveInterval(solveInterval_);
                }
                active_ = false;
                break;
            }
        }
    }
    else
    {
        Info<< "    Conditions not met" << nl;
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::write()
{
    for (auto& condition : conditions_)
    {
        condition.write();
    }

    return true;
}


// ************************************************************************* //
