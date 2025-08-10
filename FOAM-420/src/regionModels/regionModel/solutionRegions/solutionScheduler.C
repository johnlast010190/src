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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solutionRegions/solutionScheduler.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/UnbracketedTuple2/UnbracketedTuple2.H"
#include "regionProperties/regionProperties.H"

void Foam::solutionScheduler::assignUserInputValues()
{
    checkSolutionMesh();
    solutionScheduleDict_ = dict_.subDict("solutionSchedule");
    solutionMeshDict_ = dict_.subDict("solutionMesh");
    instanceNames_ = solutionScheduleDict_.toc();
    numberOfInstances_ = instanceNames_.size();
    solutionInstancePtrList_.setSize(numberOfInstances_);
    buildInstances();
}

void Foam::solutionScheduler::selectOnly
(
    const HashSet<word>& selectedInstances
)
{
    forAll(instanceNames_, inst)
    {
        if (!selectedInstances.found(instanceNames_[inst]))
        {
            solutionInstancePtrList_[inst].unselectInstance();
        }
    }
    forAll(instanceNames_, inst)
    {
        if (selectedInstances.found(instanceNames_[inst]))
        {
            currentInstanceIndex_ = inst;
            break;
        }
    }
    Info<<"Current Instance Index is "<<currentInstanceIndex_<<endl;
}

void Foam::solutionScheduler::checkSolutionMesh()
{
    if (!dict_.found("solutionMesh"))
    {
        FatalErrorInFunction
            << "Find solutionSchedule inside regionProperties but "
            << "there is no definition of solutionMesh"
            << exit(FatalError);
    }
}

void Foam::solutionScheduler::buildInstances()
{
    forAll(instanceNames_, inst)
    {
        List<List<word>> instanceRegions
        (
            solutionScheduleDict_.lookup(instanceNames_[inst])
        );
        solutionInstancePtrList_.set
        (
            inst,
            new solutionInstance
            (
                runTime_,
                instanceRegions,
                solutionMeshDict_,
                instanceNames_[inst]
            )
        );
    }
}


void Foam::solutionScheduler::checkValidInput()
{
    //TODO Write a checking input algorithm
}


void Foam::solutionScheduler::assignDefaultValues()
{
    defaultBehavior_ = true;
    instanceNames_.setSize(1);
    instanceNames_[0] = "defaultInstanceName";

    numberOfInstances_ = instanceNames_.size();
    solutionInstancePtrList_.setSize(numberOfInstances_);

    regionProperties rp(runTime_);
    label total = 0;
    forAll(rp.groupNames(), i)
    {
        total += rp[rp.groupNames()[i]].size();
    }
    labelList  regionGroup(total);
    hashedWordList regionNames;
    forAll(rp.groupNames(), i)
    {
        const wordList& regionNamesi = rp[rp.groupNames()[i]];
        forAll(regionNamesi, j)
        {
            regionNames.append(regionNamesi[j]);
        }
    }
    List<wordList> instanceRegions(1, regionNames);
    dictionary meshDict;
    forAll(regionNames, rN)
    {
        meshDict.add(regionNames[rN], regionNames[rN]);
    }
    solutionMeshDict_ = meshDict;
    solutionInstancePtrList_.set
    (
        0,
        new solutionInstance
        (
            runTime_,
            instanceRegions,
            meshDict,
            "defaultInstanceName",
            true
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionScheduler::solutionScheduler
(
    argList& args,
    const Time& runTime
)
:
    runTime_(runTime),
    dict_
    (
        IOobject
        (
            "regionProperties",
            runTime.time().constant(),
            runTime.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    regionMap_(runTime),
    meshNames_(),
    selectedInstances_(),
    numberOfInstances_(0),
    instanceCounter_(0),
    currentInstanceIndex_(0),
    solutionInstancePtrList_(0),
    defaultBehavior_(false)
{
    if
    (
        dict_.typeHeaderOk<IOdictionary>(true)
     && dict_.found("solutionSchedule")
    )
    {
        dict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        assignUserInputValues();
        if (args.optionFound("instances"))
        {
            args.optionLookup("instances")() >> selectedInstances_;
            selectOnly(selectedInstances_);
        }
        else
        {
            selectedInstances_ = instanceNames_;
        }
        Info<<"SELECTED INSTANCES "<<selectedInstances_<<endl;
    }
    else
    {
        assignDefaultValues();
    }
}


Foam::solutionScheduler::solutionScheduler
(
    const Time& runTime,
    const IOdictionary& IOdict,
    bool caseSetup
)
:
    runTime_(runTime),
    dict_
    (
        IOdict
    ),
    regionMap_(runTime),
    instanceNames_(0),
    meshNames_(),
    numberOfInstances_(0),
    instanceCounter_(0),
    currentInstanceIndex_(0),
    solutionInstancePtrList_(0),
    defaultBehavior_(false)
{
    if
    (
        dict_.typeHeaderOk<IOdictionary>(false)
     && dict_.found("solutionSchedule")
    )
    {
        dict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;

        assignUserInputValues();
    }
    else if (caseSetup)
    {
        defaultBehavior_ = true;
        instanceNames_.setSize(1);
        instanceNames_[0] = "defaultInstanceName";

        numberOfInstances_ = instanceNames_.size();
        solutionInstancePtrList_.setSize(numberOfInstances_);
    }
    else
    {
        assignDefaultValues();
    }
}

Foam::solutionScheduler::~solutionScheduler() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solutionScheduler::loop()
{
    return instanceCounter_ < numberOfInstances_;
}


const Foam::fvSolutionRegistry& Foam::solutionScheduler::getActiveSolutionRegistry
(
    const Time& runTime
)
{
    return solutionInstancePtrList_[currentInstanceIndex_].getSolutionRegistry(runTime);
}


const Foam::PtrList<Foam::fvSolutionRegistry>&
Foam::solutionScheduler::getActiveSolutionRegistryList
(
    const Time& runTime
)
{
    return solutionInstancePtrList_[currentInstanceIndex_].getSolutionRegistryList
    (
        runTime
    );
}


void Foam::solutionScheduler::operator++()
{
    bool foundNextInstance= false;
    for (int i=0; i<numberOfInstances_; ++i)
    {
        nextInstanceIndex_ = (currentInstanceIndex_+1+i)%numberOfInstances_;
        instanceCounter_++;
        if (selectedInstances_.found(instanceNames_[nextInstanceIndex_]))
        {
            foundNextInstance = true;
            break;
        }
    }

    if (foundNextInstance)
    {
        solutionInstance& currentInstance = solutionInstancePtrList_[currentInstanceIndex_];
        const solutionInstance& nextInstance = solutionInstancePtrList_[nextInstanceIndex_];
        const wordList& currentMeshNames = currentInstance.meshNames();
        const wordList& nextMeshNames = nextInstance.meshNames();
        forAll(currentMeshNames, meshi)
        {
            bool deleteMesh = true;
            forAll(nextMeshNames, nextMeshi)
            {
                if (currentMeshNames[meshi] == nextMeshNames[nextMeshi])
                {
                    deleteMesh = false;
                }
            }
            if (deleteMesh)
            {
                currentInstance.releaseUnusedMeshes
                (
                    runTime_,
                    currentMeshNames[meshi]
                );
            }
            currentInstance.deactivate();
        }
        currentInstanceIndex_ = nextInstanceIndex_;
    }
}


void Foam::solutionScheduler::operator++(int)
{
    return operator++();
}


Foam::dictionary Foam::solutionScheduler::getSchedulerDictionary() const
{
    dictionary dict;
    if (defaultBehavior_)
    {
        return dict;
    }
    return solutionScheduleDict_;
}


Foam::dictionary Foam::solutionScheduler::getSolutionMeshDictionary() const
{
    return solutionMeshDict_;
}
