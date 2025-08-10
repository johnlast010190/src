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
    (c) 2010-2020 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/Time/Time.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include "global/argList/argList.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "containers/HashTables/HashSet/HashSet.H"

#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

#include "include/demandDrivenData.H"


#include <sstream>

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(loadLibs, 0);
    defineTypeNameAndDebug(Time, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::Time::stopAtControls,
        5
    >::names[] =
    {
        "endTime",
        "noWriteNow",
        "writeNow",
        "nextWrite",
        "clockTime"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::Time::writeControls,
        7
    >::names[] =
    {
        "timeStep",
        "runTime",
        "adjustableRunTime",
        "clockTime",
        "cpuTime",
        "baselineTime",
        "optimisationCycle"
    };
}

const Foam::NamedEnum<Foam::Time::stopAtControls, 5>
    Foam::Time::stopAtControlNames_;

const Foam::NamedEnum<Foam::Time::writeControls, 7>
    Foam::Time::writeControlNames_;

Foam::Time::fmtflags Foam::Time::format_(Foam::Time::general);

int Foam::Time::precision_(6);

const int Foam::Time::maxPrecision_(3 - log10(SMALL));

Foam::word Foam::Time::controlDictName("controlDict");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Time::adjustDeltaT()
{
    bool adjustTime = false;
    scalar timeToNextWrite = VGREAT;

    if (writeControl_ == wcAdjustableRunTime)
    {
        adjustTime = true;
        timeToNextWrite = max
        (
            0.0,
            (writeTimeIndex_ + 1)*writeInterval_ - (value() - startTime_)
        );
    }

    if (adjustTime)
    {
        scalar nSteps = timeToNextWrite/deltaT_;

        // For tiny deltaT the label can overflow!
        if (nSteps < labelMax)
        {
            // nSteps can be < 1 so make sure at least 1
            label nStepsToNextWrite = max(1, round(nSteps));

            scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

            // Control the increase of the time step to within a factor of 2
            // and the decrease within a factor of 5.
            if (newDeltaT >= deltaT_)
            {
                deltaT_ = min(newDeltaT, 2.0*deltaT_);
            }
            else
            {
                deltaT_ = max(newDeltaT, 0.2*deltaT_);
            }
        }
    }

    functionObjects_.adjustTimeStep();
}


void Foam::Time::setControls()
{
    // Make sure we have the correct time-format parameters
    controlDict_.readIfPresent("timePrecision", precision_);

    // default is to resume calculation from "latestTime"
    const word startFrom = controlDict_.lookupOrDefault<word>
    (
        "startFrom",
        "latestTime"
    );

    if (startFrom == "startTime")
    {
        controlDict_.lookup("startTime") >> startTime_;
    }
    else
    {
        // Search directory for valid time directories
        instantList timeDirs = findTimes(path(), constant());

        if (startFrom == "firstTime")
        {
            if (timeDirs.size())
            {
                if (timeDirs[0].name() == constant() && timeDirs.size() >= 2)
                {
                    startTime_ = timeDirs[1].value();
                }
                else
                {
                    startTime_ = timeDirs[0].value();
                }
            }
        }
        else if (startFrom == "latestTime")
        {
            if (timeDirs.size())
            {
                startTime_ = timeDirs.last().value();
            }
        }
        else
        {
            FatalIOErrorInFunction(controlDict_)
                << "expected startTime, firstTime or latestTime"
                << " found '" << startFrom << "'"
                << exit(FatalIOError);
        }
    }

    // Set time to get correct timePath()
    setTime(startTime_, 0);

    readDict();

    // Check if time directory exists
    // If not increase time precision to see if it is formatted differently.
    if (!fileHandler().exists(timePath(), false))
    {
        int oldPrecision = precision_;
        int requiredPrecision = -1;
        bool found = false;
        word oldTime(timeName());

        for
        (
            precision_ = maxPrecision_;
            precision_ > oldPrecision;
            precision_--
        )
        {
            // Update the time formatting
            setTime(startTime_, 0);

            // Check the existence of the time directory with the new format
            word newTime(timeName());

            // If already found and the time name changed, don't go further
            // If it didn't change, keep going to find the smallest precision
            if (found)
            {
                if (newTime == oldTime)
                {
                    requiredPrecision = precision_;
                }
                else
                {
                    break;
                }
            }
            else if (fileHandler().exists(timePath(), false))
            {
                requiredPrecision = precision_;
                found = true;
            }

            oldTime = newTime;
        }

        if (requiredPrecision > 0)
        {
            // Update the time precision
            precision_ = requiredPrecision;

            // Update the time formatting
            setTime(startTime_, 0);

            WarningInFunction
                << "Increasing the timePrecision from " << oldPrecision
                << " to " << precision_
                << " to support the formatting of the current time directory "
                << timeName() << nl << endl;
        }
        else
        {
            // Could not find time directory so assume it is not present
            precision_ = oldPrecision;

            // Revert the time formatting
            setTime(startTime_, 0);
        }
    }

    if (Pstream::parRun())
    {
        scalar sumStartTime = startTime_;
        reduce(sumStartTime, sumOp<scalar>());
        if
        (
            mag(Pstream::nProcs()*startTime_ - sumStartTime)
          > Pstream::nProcs()*deltaT_/10.0
        )
        {
            FatalIOErrorInFunction(controlDict_)
                << "Start time is not the same for all processors" << nl
                << "processor " << Pstream::myProcNo() << " has startTime "
                << startTime_ << exit(FatalIOError);
        }
    }

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (timeDict.readIfPresent("index", startTimeIndex_))
    {
        timeIndex_ = startTimeIndex_;
    }
    setTime(startTime_, startTimeIndex_);

    // Re-initialise table deltaT now that timeIndex_ is set
    if (deltaTDataPtr_.valid())
    {
        if (deltaTcontrol_ == "index")
        {
            deltaT_ = deltaTDataPtr_->value(timeIndex_);
        }
        else if (deltaTcontrol_ == "time")
        {
            deltaT_ = deltaTDataPtr_->value(value());
        }
        deltaTSave_ = deltaT_;
        deltaT0_ = deltaT_;
    }

    // Read and set the deltaT only if time-step adjustment is active
    // otherwise use the deltaT from the controlDict
    if (controlDict_.lookupOrDefault<Switch>("adjustTimeStep", false))
    {
        if (timeDict.readIfPresent("deltaT", deltaT_))
        {
            deltaTSave_ = deltaT_;
            deltaT0_ = deltaT_;
        }
    }

    timeDict.readIfPresent("deltaT0", deltaT0_);


    // Check if values stored in time dictionary are consistent

    // 1. Based on time name
    bool checkValue = true;

    string storedTimeName;
    if (timeDict.readIfPresent("name", storedTimeName))
    {
        if (storedTimeName == timeName())
        {
            // Same time. No need to check stored value
            checkValue = false;
        }
    }

    // 2. Based on time value
    //    (consistent up to the current time writing precision so it won't
    //     trigger if we just change the write precision)
    if (checkValue)
    {
        scalar storedTimeValue;
        if (timeDict.readIfPresent("value", storedTimeValue))
        {
            word storedTimeName(timeName(storedTimeValue));

            if (storedTimeName != timeName())
            {
                IOWarningInFunction(timeDict)
                    << "Time read from time dictionary " << storedTimeName
                    << " differs from actual time " << timeName() << '.' << nl
                    << "    This may cause unexpected database behaviour."
                    << " If you are not interested" << nl
                    << "    in preserving time state delete"
                    << " the time dictionary."
                    << endl;
            }
        }
    }
}


void Foam::Time::setMonitoring(const bool forceProfiling)
{

#if !defined( WIN32 ) && !defined( WIN64 )
    const dictionary* profilingDict = controlDict_.subDictPtr("profiling");
    if (!profilingDict)
    {
        // ... or from etc/controlDict
        profilingDict = debug::controlDict().subDictPtr("profiling");
    }

    // initialize profiling on request
    // otherwise rely on profiling entry within controlDict
    // and skip if 'active' keyword is explicitly set to false
    if (forceProfiling)
    {
        profiling::initialize
        (
            IOobject
            (
                "profiling",
                timeName(),
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
    else if
    (
        profilingDict
     && profilingDict->lookupOrDefault<bool>("active", true)
    )
    {
        profiling::initialize
        (
            *profilingDict,
            IOobject
            (
                "profiling",
                timeName(),
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
#endif

    // Time objects not registered so do like objectRegistry::checkIn ourselves.
    if (runTimeModifiable_)
    {
        // Monitor all files that controlDict depends on
        fileHandler().addWatches(controlDict_, controlDict_.files());
    }

    // Clear dependent files - not needed now
    controlDict_.files().clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Time::Time
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

#if !defined( WIN32 ) && !defined( WIN64 )
    loopProfiling_(nullptr),
#endif

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),
    endClockTime_(0),
    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(2E9),
    writeEndTime_(true),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),
    baselineTime_(0),

    functionObjects_(*this, enableFunctionObjects)
{
    libs().open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();
    setMonitoring();
}


Foam::Time::Time
(
    const word& controlDictName,
    const argList& args,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        args.parRunControl().parRun(),
        args.rootPath(),
        args.distributed(),
        args.globalCaseName(),
        args.caseName(),
        systemName,
        constantName
    ),

    objectRegistry(*this),

#if !defined( WIN32 ) && !defined( WIN64 )
    loopProfiling_(nullptr),
#endif

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),
    endClockTime_(0),
    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(2E9),
    writeEndTime_(true),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),
    baselineTime_(0),

    functionObjects_
    (
        *this,
        argList::validOptions.found("withFunctionObjects")
      ? args.optionFound("withFunctionObjects")
      : argList::validOptions.found("noFunctionObjects")
      ? !args.optionFound("noFunctionObjects")
      : false
    )
{
    libs().open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();

    // '-profiling' = force profiling, ignore controlDict entry
    setMonitoring(args.optionFound("profiling"));

}


Foam::Time::Time
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

#if !defined( WIN32 ) && !defined( WIN64 )
    loopProfiling_(nullptr),
#endif

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dict
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),
    endClockTime_(0),
    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(2E9),
    writeEndTime_(true),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),
    baselineTime_(0),

    functionObjects_(*this, enableFunctionObjects)
{
    libs().open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    // Since could not construct regIOobject with setting:
    controlDict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();
    setMonitoring();
}


Foam::Time::Time
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

#if !defined( WIN32 ) && !defined( WIN64 )
    loopProfiling_(nullptr),
#endif

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),
    endClockTime_(0),
    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(2E9),
    writeEndTime_(true),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),
    baselineTime_(0),

    functionObjects_(*this, enableFunctionObjects)
{
    libs().open(controlDict_, "libs");
    setMonitoring(); // for profiling etc
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Time::~Time()
{
#if !defined( WIN32 ) && !defined( WIN64 )
    deleteDemandDrivenData(loopProfiling_);
#endif

    forAllReverse(controlDict_.watchIndices(), i)
    {
        fileHandler().removeWatch(controlDict_.watchIndices()[i]);
    }

    // Destroy function objects first
    functionObjects_.clear();

#if !defined( WIN32 ) && !defined( WIN64 )
    // Clean up profiling
    profiling::stop(*this);
#endif

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::Time::timeName(const scalar t, const int precision)
{
    std::ostringstream buf;
    buf.setf(ios_base::fmtflags(format_), ios_base::floatfield);
    buf.precision(precision);
    buf << t;
    return buf.str();
}


Foam::word Foam::Time::timeName() const
{
    return dimensionedScalar::name();
}


Foam::instantList Foam::Time::times() const
{
    return findTimes(path(), constant());
}


Foam::word Foam::Time::findInstancePath
(
    const fileName& directory,
    const instant& t
) const
{
    // Simplified version: use findTimes (readDir + sort). The expensive
    // bit is the readDir, not the sorting. Tbd: avoid calling findInstancePath
    // from filePath.

    instantList timeDirs = findTimes(path(), constant());
    // Note:
    // - times will include constant (with value 0) as first element.
    //   For backwards compatibility make sure to find 0 in preference
    //   to constant.
    // - list is sorted so could use binary search

    forAllReverse(timeDirs, i)
    {
        if (t.equal(timeDirs[i].value()))
        {
            return timeDirs[i].name();
        }
    }

    return word::null;
}


Foam::word Foam::Time::findInstancePath(const instant& t) const
{
    return findInstancePath(path(), t);
}


Foam::instant Foam::Time::findClosestTime(const scalar t) const
{
    instantList timeDirs = findTimes(path(), constant());

    // There is only one time (likely "constant") so return it
    if (timeDirs.size() == 1)
    {
        return timeDirs[0];
    }

    if (t < timeDirs[1].value())
    {
        return timeDirs[1];
    }
    else if (t > timeDirs.last().value())
    {
        return timeDirs.last();
    }

    label nearestIndex = -1;
    scalar deltaT = GREAT;

    for (label timei=1; timei < timeDirs.size(); ++timei)
    {
        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return timeDirs[nearestIndex];
}


Foam::label Foam::Time::findClosestTimeIndex
(
    const instantList& timeDirs,
    const scalar t,
    const word& constantName
)
{
    label nearestIndex = -1;
    scalar deltaT = GREAT;

    forAll(timeDirs, timei)
    {
        if (timeDirs[timei].name() == constantName) continue;

        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return nearestIndex;
}


Foam::label Foam::Time::restartTimeIndex() const
{
    return (timeIndex() - startTimeIndex());
}

Foam::label Foam::Time::startTimeIndex() const
{
    return startTimeIndex_;
}


Foam::dimensionedScalar Foam::Time::startTime() const
{
    return dimensionedScalar("startTime", dimTime, startTime_);
}


Foam::dimensionedScalar Foam::Time::endTime() const
{
    return dimensionedScalar("endTime", dimTime, endTime_);
}


bool Foam::Time::run() const
{
#if !defined( WIN32 ) && !defined( WIN64 )
    deleteDemandDrivenData(loopProfiling_);
#endif

    bool isRunning = value() < (endTime_ - 0.5*deltaT_);

    if (!subCycling_)
    {
        // Only execute when the condition is no longer true
        // ie, when exiting the control loop
        if (!isRunning && timeIndex_ != startTimeIndex_)
        {
            // Ensure functionObjects execute on last time step
            // (and hence write uptodate functionObjectProperties)
#if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(foExec, "functionObjects.execute()");
#endif

            functionObjects_.execute();

#if !defined( WIN32 ) && !defined( WIN64 )
            endProfiling(foExec);
            addProfiling(foEnd, "functionObjects.end()");
#endif
            functionObjects_.end();

#if !defined( WIN32 ) && !defined( WIN64 )
            endProfiling(foEnd);
#endif
        }
    }

    if (isRunning)
    {
        if (!subCycling_)
        {
            const_cast<Time&>(*this).readModifiedObjects();

            if (timeIndex_ == startTimeIndex_)
            {
#if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(functionObjects, "functionObjects.start()");
#endif
                functionObjects_.start();
            }
            else
            {
#if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(functionObjects, "functionObjects.execute()");
#endif
                functionObjects_.execute();
            }

            // Check if the execution of functionObjects require re-reading
            // any files. This moves effect of e.g. 'timeActivatedFileUpdate'
            // one time step forward. Note that we cannot call
            // readModifiedObjects from within timeActivatedFileUpdate since
            // it might re-read the functionObjects themselves (and delete
            // the timeActivatedFileUpdate one)
            if (functionObjects_.filesModified())
            {
                const_cast<Time&>(*this).readModifiedObjects();
            }
        }

        // Update the "is-running" status following the
        // possible side-effects from functionObjects
        isRunning = value() < (endTime_ - 0.5*deltaT_);

#if !defined( WIN32 ) && !defined( WIN64 )
        // (re)trigger profiling
        if (profiling::active())
        {
            loopProfiling_ =
                new profilingTrigger("time.run() " + objectRegistry::name());
        }
#endif
    }

    return isRunning;
}


bool Foam::Time::loop()
{
    const bool isRunning = run();

    if (isRunning)
    {
        operator++();
    }

    return isRunning;
}


bool Foam::Time::end() const
{
    return value() > (endTime_ + 0.5*deltaT_);
}


bool Foam::Time::stopAt(const stopAtControls sa) const
{
    const bool changed = (stopAt_ != sa);
    stopAt_ = sa;

    // adjust endTime
    if (sa == saEndTime)
    {
        controlDict_.lookup("endTime") >> endTime_;
    }
    else if (sa == saClockTime)
    {
        controlDict_.lookup("endTime") >> endTime_;
        controlDict_.lookup("endClockTime") >> endClockTime_;
    }
    else
    {
        endTime_ = GREAT;
    }
    return changed;
}


void Foam::Time::setTime(const Time& t)
{
    value() = t.value();
    dimensionedScalar::name() = t.dimensionedScalar::name();
    timeIndex_ = t.timeIndex_;
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const instant& inst, const label newIndex)
{
    value() = inst.value();
    dimensionedScalar::name() = inst.name();
    timeIndex_ = newIndex;

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    timeDict.readIfPresent("deltaT", deltaT_);
    timeDict.readIfPresent("deltaT0", deltaT0_);
    timeDict.readIfPresent("index", timeIndex_);
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const dimensionedScalar& newTime, const label newIndex)
{
    setTime(newTime.value(), newIndex);
}


void Foam::Time::setTime(const scalar newTime, const label newIndex)
{
    value() = newTime;
    dimensionedScalar::name() = timeName(timeToUserTime(newTime));
    timeIndex_ = newIndex;
    fileHandler().setTime(*this);
}


void Foam::Time::setEndTime(const dimensionedScalar& endTime)
{
    setEndTime(endTime.value());
}


void Foam::Time::setEndTime(const scalar endTime)
{
    endTime_ = endTime;
}


void Foam::Time::setDeltaT
(
    const dimensionedScalar& deltaT,
    const bool adjust
)
{
    setDeltaT(deltaT.value(), adjust);
}


void Foam::Time::setDeltaT(const scalar deltaT, const bool adjust)
{
    if (deltaT_ != deltaT)
    {
        deltaT_ = deltaT;
        deltaTchanged_ = true;

        if (adjust)
        {
            adjustDeltaT();
        }
    }
}


Foam::TimeState Foam::Time::subCycle(const label nSubCycles)
{
    subCycling_ = true;
    prevTimeState_.set(new TimeState(*this));

    setTime(*this - deltaT(), (timeIndex() - 1)*nSubCycles);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;
    deltaTSave_ = deltaT0_;

    return prevTimeState();
}


void Foam::Time::endSubCycle()
{
    if (subCycling_)
    {
        TimeState::operator=(prevTimeState());
        prevTimeState_.clear();
        subCycling_ = false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Time& Foam::Time::operator+=(const dimensionedScalar& deltaT)
{
    return operator+=(deltaT.value());
}


Foam::Time& Foam::Time::operator+=(const scalar deltaT)
{
    setDeltaT(deltaT);
    return operator++();
}


Foam::Time& Foam::Time::operator++()
{
    if (deltaTDataPtr_.valid() && !subCycling_)
    {
        scalar ddt = deltaT_;
        if (deltaTcontrol_ == "index")
        {
            setDeltaT(deltaTDataPtr_->value(timeIndex_ + 1));
        }
        else if (deltaTcontrol_ == "time")
        {
            setDeltaT(deltaTDataPtr_->value(value() + deltaT_));
        }
        ddt -= deltaT_;
        ddt *= -1;

        Info<< "Function1 type: " << deltaTDataPtr_->type()
             << " - iter: " << timeIndex_ + 1
             << " - ddt: " << ddt
             << " - deltaT: " << deltaT_
             << endl;
    }

    deltaT0_ = deltaTSave_;
    deltaTSave_ = deltaT_;

    // Save old time value and name
    const scalar oldTimeValue = timeToUserTime(value());
    const word oldTimeName = dimensionedScalar::name();

    // Increment time
    setTime(value() + deltaT_, timeIndex_ + 1);

    if (!subCycling_)
    {
        // If the time is very close to zero reset to zero
        if (mag(value()) < 10*SMALL*deltaT_)
        {
            setTime(0.0, timeIndex_);
        }

        if (stopAt_ == saClockTime)
        {
            label currClockTime = returnReduce
            (
                label(elapsedClockTime()), maxOp<label>()
            );
            if (currClockTime >= endClockTime_)
            {
                stopAt_ = saWriteNow;
            }
        }

        if (sigStopAtWriteNow_.active() || sigWriteNow_.active())
        {
            // A signal might have been sent on one processor only
            // Reduce so all decide the same.

            label flag = 0;
            if (sigStopAtWriteNow_.active() && stopAt_ == saWriteNow)
            {
                flag += 1;
            }
            if (sigWriteNow_.active() && writeOnce_)
            {
                flag += 2;
            }
            reduce(flag, maxOp<label>());

            if (flag & 1)
            {
                stopAt_ = saWriteNow;
            }
            if (flag & 2)
            {
                writeOnce_ = true;
            }
        }

        writeTime_ = false;

        switch (writeControl_)
        {
            case wcTimeStep:
                writeTime_ = !(timeIndex_ % label(writeInterval_));
            break;

            case wcRunTime:
            case wcAdjustableRunTime:
            {
                label writeIndex = label
                (
                    ((value() - startTime_) + 0.5*deltaT_)
                  / writeInterval_
                );

                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcCpuTime:
            {
                label writeIndex = label
                (
                    returnReduce(elapsedCpuTime(), maxOp<double>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcClockTime:
            {
                label writeIndex = label
                (
                    returnReduce(label(elapsedClockTime()), maxOp<label>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcBaselineTime:
            {
                scalar tStart = value() - deltaT_;

                if (value() >= baselineTime_ && tStart >= baselineTime_)
                {
                    label outputIndexStart = label
                    (
                        ((tStart - baselineTime_))
                        / writeInterval_
                    );

                    label outputIndexEnd = label
                    (
                        (value() - baselineTime_)
                        / writeInterval_
                    );

                    if (outputIndexEnd > outputIndexStart)
                    {
                        writeTime_ = true;
                        writeTimeIndex_ = outputIndexEnd;
                    }
                }
            }
            break;

            case wcOptimizationCycle:
            {
                int intValue = round (value());
                const scalar timeTol =
                    max(min(pow(10.0, -precision_), 0.1*deltaT_), SMALL);
                if (mag(value()- intValue) < timeTol)
                {
                    writeTime_ = !(intValue % int (writeInterval_));
                }
            }
            break;
        }


        // Check if endTime needs adjustment to stop at the next run()/end()
        if (!end())
        {
            if (stopAt_ == saNoWriteNow)
            {
                endTime_ = value();
            }
            else if (stopAt_ == saWriteNow)
            {
                endTime_ = value();
                writeTime_ = true;
            }
            else if (stopAt_ == saNextWrite && writeTime_ == true)
            {
                endTime_ = value();
            }
        }

        if (writeEndTime_ && value() > (endTime_ - 0.5*deltaT_))
        {
            //write fields if final time-step and end time write enabled
            writeTime_ = true;
        }

        // Override writeTime if one-shot writing
        if (writeOnce_)
        {
            writeTime_ = true;
            writeOnce_ = false;
        }

        // Adjust the precision of the time directory name if necessary
        if (writeTime_)
        {
            // Tolerance used when testing time equivalence
            const scalar timeTol =
                max(min(pow(10.0, -precision_), 0.1*deltaT_), SMALL);

            // User-time equivalent of deltaT
            const scalar userDeltaT = timeToUserTime(deltaT_);

            // Time value obtained by reading timeName
            scalar timeNameValue = -VGREAT;

            // Check that new time representation differs from old one
            // reinterpretation of the word
            if
            (
                readScalar(dimensionedScalar::name().c_str(), timeNameValue)
             && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
            )
            {
                int oldPrecision = precision_;
                while
                (
                    precision_ < maxPrecision_
                 && readScalar(dimensionedScalar::name().c_str(), timeNameValue)
                 && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
                )
                {
                    precision_++;
                    setTime(value(), timeIndex());
                }

                if (precision_ != oldPrecision)
                {
                    WarningInFunction
                        << "Increased the timePrecision from " << oldPrecision
                        << " to " << precision_
                        << " to distinguish between timeNames at time "
                        << dimensionedScalar::name()
                        << endl;

                    if (precision_ == maxPrecision_)
                    {
                        // Reached maxPrecision limit
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << nl
                            << "    The maximum time precision has been reached"
                               " which might result in overwriting previous"
                               " results."
                            << endl;
                    }

                    // Check if round-off error caused time-reversal
                    scalar oldTimeNameValue = -VGREAT;
                    if
                    (
                        readScalar(oldTimeName.c_str(), oldTimeNameValue)
                     && (
                            sign(timeNameValue - oldTimeNameValue)
                         != sign(deltaT_)
                        )
                    )
                    {
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << " is set to an instance prior to the "
                               "previous one "
                            << oldTimeName << nl
                            << "    This might result in temporal "
                               "discontinuities."
                            << endl;
                    }
                }
            }
        }
    }

    return *this;
}


Foam::Time& Foam::Time::operator++(int)
{
    return operator++();
}

void Foam::Time::setIOStream
(
    const IOstream::streamFormat    fmt,
    const IOstream::versionNumber   ver,
    const IOstream::compressionType cmp
)
{
    writeFormat_ = fmt;
    writeVersion_ = ver;
    writeCompression_ = cmp;
}
// ************************************************************************* //
