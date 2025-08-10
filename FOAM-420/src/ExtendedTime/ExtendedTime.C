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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "ExtendedTime.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ExtendedTime, 0);

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

bool
ExtendedTime::inverseRun() const
{
    bool running = value() > (initialTime_ + 0.5*primalDeltaT_);

    functionObjectList& foList =
        const_cast<functionObjectList&>(functionObjects());

    if (!subCycling_)
    {
        const_cast<ExtendedTime&>(*this).readModifiedObjects();

        forAll(functionObjectsInReverse_, foI)
        {
            label foID = foList.findObjectID(functionObjectsInReverse_[foI]);

            if
            (
                foID != -1
             && value() < finalTime_ - 0.5*reverseDeltaT_
            )
            {
                foList[foID].execute();
            }
        }
    }

    // Update the "running" status following the
    // possible side-effects from functionObjects
    running = value() > initialTime_ + 0.5*primalDeltaT_;

    return running;
}

bool
ExtendedTime::inverseLoop()
{
    bool running = inverseRun();

    if (running)
    {
        operator--();
    }

    return running;
}

bool
ExtendedTime::inverseEnd() const
{
    return value() > initialTime_ + 0.5*primalDeltaT_;
}

void
ExtendedTime::init()
{
    cpTable_ =
        readList<word>
        (
            IStringStream
            (
               "(U p phi nut nuSgs nuTilda k omega)"
            )()
        );
    cpTable_.append
    (
        readList<word>
        (
            IStringStream
            (
               "(U_0 phi_0 nuTilda_0 k_0 omega_0)"
            )()
        )
    );

    IOdictionary adjointProperties
    (
        IOobject
        (
            "adjointProperties",
            caseSystem(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    nCP_ = adjointProperties.subDict("solution").
        lookupOrDefault<label>
        (
            "nCheckpoints",
            1
        );

    functionObjectsInReverse_ = adjointProperties.subDict("solution").
        lookupOrDefault<wordList>
        (
            "FOsInReverse",
            wordList(0)
        );

    reverseTimestep_ = adjointProperties.subDict("solution").
        lookupOrDefault<label>
        (
            "reverseTimestep",
            1
        );

    primalDeltaT_ = deltaT_;
    reverseDeltaT_ = reverseTimestep_*deltaT_;

    if (adjointProperties.subDict("solution").found("startTime"))
    {
        scalar initTime =
            readScalar
            (
                adjointProperties.subDict("solution").
                   lookup("startTime")
            );

        if (initTime > startTime_ - 0.5*deltaT_)
            initialTime_ = initTime;
        else
            FatalErrorInFunction
                << "startTime for adjoint smaller than the job's startTime"
                << exit(FatalError);
    }

    checkpoint_.resize(nCP_);
    cpTimeValue_.resize(nCP_);
    cpDeltaTValue_.resize(nCP_);
    cpDeltaT0Value_.resize(nCP_);
    cpTimeIndex_.resize(nCP_);
    cpLevel_.resize(nCP_);
    cpAddressing_.resize(nCP_);

    storeCPs_ = adjointProperties.subDict("solution").
        lookupOrDefault<Switch>
        (
            "storeCPsAtEnd",
            true
        );

    forAll(checkpoint_, cpi)
    {
        checkpoint_.set
        (
            cpi,
            new objectRegistry
            (
                IOobject
                (
                    "checkpoint_[" + word(std::to_string(cpi)) + "]",
                    timeName(),
                    "../checkpoints",
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        cpLevel_[cpi] = -1;
        cpAddressing_[cpi] = cpi;
    }
}

void
ExtendedTime::addCP()
{
    if
    (
        value() + reverseDeltaT_ > (initialTime_ - 0.5*primalDeltaT_)
     && value() < (endTime_ - 0.5*primalDeltaT_)
     && isValidTimeStepInReverse(value())
    )
    {
        label numCPsNow = 0;
        forAll(checkpoint_, cpi)
        {
            if (cpLevel_[cpAddressing_[cpi]] > -1)
                numCPsNow++;
        }

        label idis = 0;
        label maxLevel = 0;
        for (int i = nCP_ - 1; i >= 1; i--)
        {
            if (maxLevel > cpLevel_[cpAddressing_[i]])
            {
                if (idis == 0)
                    idis = i;
            }
            else
                maxLevel = cpLevel_[cpAddressing_[i]];
        }

        if (numCPsNow < nCP_)
        {
            storeCP(cpAddressing_[numCPsNow]);
            if (numCPsNow > 1)
                cpLevel_[cpAddressing_[numCPsNow]] = 0;
            else if (numCPsNow == 1)
                cpLevel_[cpAddressing_[numCPsNow]] = 999;
            else
                cpLevel_[cpAddressing_[numCPsNow]] = 999;
        }
        else if (idis > 0)
        {
            label indirect_pnt = cpAddressing_[idis];

            for (label i = idis; i < nCP_ - 1; i++)
            {
                cpAddressing_[i] = cpAddressing_[i + 1];
            }
            cpAddressing_[nCP_ - 1] = indirect_pnt;

            storeCP(cpAddressing_[nCP_ - 1]);
            cpLevel_[cpAddressing_[nCP_ - 1]] = 0;
        }
        else
        {
            storeCP(cpAddressing_[nCP_ - 1]);
            cpLevel_[cpAddressing_[nCP_ - 1]] += 1;
        }
    }
}

void
ExtendedTime::storeCP(label cpi)
{
    addCurrentFieldsToCP<scalar, fvPatchField, volMesh>(cpi);
    addCurrentFieldsToCP<vector, fvPatchField, volMesh>(cpi);
    addCurrentFieldsToCP<tensor, fvPatchField, volMesh>(cpi);
    addCurrentFieldsToCP<symmTensor, fvPatchField, volMesh>(cpi);
    addCurrentFieldsToCP<sphericalTensor, fvPatchField, volMesh>(cpi);
    addCurrentFieldsToCP<scalar, fvsPatchField, surfaceMesh>(cpi);
    addCurrentFieldsToCP<vector, fvsPatchField, surfaceMesh>(cpi);
    addCurrentFieldsToCP<tensor, fvsPatchField, surfaceMesh>(cpi);
    addCurrentFieldsToCP<symmTensor, fvsPatchField, surfaceMesh>(cpi);
    addCurrentFieldsToCP<sphericalTensor, fvsPatchField, surfaceMesh>(cpi);

    cpTimeIndex_[cpi] = timeIndex();
    cpTimeValue_[cpi] = value();
    cpDeltaTValue_[cpi] = deltaT_;
    cpDeltaT0Value_[cpi] = deltaT0_;
}

bool
ExtendedTime::cpIsAvailable()
{
    bool found = false;

    forAll(checkpoint_, cpi)
    {
        label cp = cpAddressing_[cpi];
        if
        (
            cpLevel_[cp] > -1
         && cpTimeIndex_[cp] == timeIndex()
        )
        {
            found = true;
        }
        else if (cpTimeIndex_[cp] > timeIndex())
            cpLevel_[cp] = -1;
    }

    return found;
}

void
ExtendedTime::restoreCurrentTimeCP()
{
    bool found = false;

    forAll(checkpoint_, cpi)
    {
        label cp = cpAddressing_[cpi];
        if
        (
            cpLevel_[cp] > -1
         && cpTimeIndex_[cp] == timeIndex()
         && !found
        )
        {
            restoreCP(cp);
            cpLevel_[cp] = -1;
            found = true;
        }
    }
    if (!found)
    {
//TODO        Info<< "Current primal checkpoint not found" << endl;
        FatalError.exit();
    }
}

void
ExtendedTime::restoreLatestTimeCP()
{
    label numCPsNow = 0;

    forAll(checkpoint_, cpi)
    {
        if (cpLevel_[cpi] > -1)
            numCPsNow++;
    }

    label lastCP = 0;
    if (numCPsNow != 0)
        lastCP  = cpAddressing_[numCPsNow - 1];

    if (cpTimeIndex_[lastCP] > timeIndex())
    {
        cpLevel_[lastCP] = -1;
        numCPsNow--;
        lastCP = cpAddressing_[numCPsNow - 1];
    }

    restoreCP(lastCP);
}

void
ExtendedTime::restoreCP(label cpi)
{
    setTime(cpTimeValue_[cpi], cpTimeIndex_[cpi]);
    deltaT_ = cpDeltaTValue_[cpi];
    deltaT0_ = cpDeltaT0Value_[cpi];

    restoreFieldsFromCP<scalar, fvPatchField, volMesh>(cpi);
    restoreFieldsFromCP<vector, fvPatchField, volMesh>(cpi);
    restoreFieldsFromCP<tensor, fvPatchField, volMesh>(cpi);
    restoreFieldsFromCP<symmTensor, fvPatchField, volMesh>(cpi);
    restoreFieldsFromCP<sphericalTensor, fvPatchField, volMesh>(cpi);
    restoreFieldsFromCP<scalar, fvsPatchField, surfaceMesh>(cpi);
    restoreFieldsFromCP<vector, fvsPatchField, surfaceMesh>(cpi);
    restoreFieldsFromCP<tensor, fvsPatchField, surfaceMesh>(cpi);
    restoreFieldsFromCP<symmTensor, fvsPatchField, surfaceMesh>(cpi);
    restoreFieldsFromCP<sphericalTensor, fvsPatchField, surfaceMesh>(cpi);
}


bool
ExtendedTime::isValidTimeStepInReverse(const scalar timeValue)
{
    return
        (
            fmod(finalTime_ - timeValue + 0.25*primalDeltaT_, reverseDeltaT_)
          < 0.5*primalDeltaT_
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ExtendedTime::ExtendedTime
(
    const word& name,
    const argList& args,
    const word& systemName,
    const word& constantName
)
:
    Time(name, args),
    initialTime_(startTime_),
    finalTime_(endTime_),
    functionObjectsInReverse_(0),
    reverseTimestep_(1),
    primalDeltaT_(0.0),
    reverseDeltaT_(0.0),
    reverseDeltaTSave_(0.0),
    cpTable_(0),
    nCP_(0),
    storeCPs_(true),
    checkpoint_(nCP_),
    cpTimeValue_(nCP_),
    cpDeltaTValue_(nCP_),
    cpDeltaT0Value_(nCP_),
    cpTimeIndex_(nCP_),
    cpLevel_(nCP_),
    cpAddressing_(nCP_)
{
    init();
}

ExtendedTime::ExtendedTime
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time(name, rootPath, caseName),
    initialTime_(startTime_),
    finalTime_(endTime_),
    functionObjectsInReverse_(0),
    reverseTimestep_(1),
    primalDeltaT_(0.0),
    reverseDeltaT_(0.0),
    reverseDeltaTSave_(0.0),
    cpTable_(0),
    nCP_(0),
    storeCPs_(true),
    checkpoint_(nCP_),
    cpTimeValue_(nCP_),
    cpDeltaTValue_(nCP_),
    cpDeltaT0Value_(nCP_),
    cpTimeIndex_(nCP_),
    cpLevel_(nCP_),
    cpAddressing_(nCP_)
{
    init();
}

ExtendedTime::ExtendedTime
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time(dict, rootPath, caseName),
    initialTime_(startTime_),
    finalTime_(endTime_),
    functionObjectsInReverse_(0),
    reverseTimestep_(1),
    primalDeltaT_(0.0),
    reverseDeltaT_(0.0),
    reverseDeltaTSave_(0.0),
    cpTable_(0),
    nCP_(0),
    storeCPs_(true),
    checkpoint_(nCP_),
    cpTimeValue_(nCP_),
    cpDeltaTValue_(nCP_),
    cpDeltaT0Value_(nCP_),
    cpTimeIndex_(nCP_),
    cpLevel_(nCP_),
    cpAddressing_(nCP_)
{
    init();
}

ExtendedTime::ExtendedTime
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time(rootPath, caseName),
    initialTime_(startTime_),
    finalTime_(endTime_),
    functionObjectsInReverse_(0),
    reverseTimestep_(1),
    primalDeltaT_(0.0),
    reverseDeltaT_(0.0),
    reverseDeltaTSave_(0.0),
    cpTable_(0),
    nCP_(0),
    storeCPs_(true),
    checkpoint_(nCP_),
    cpTimeValue_(nCP_),
    cpDeltaTValue_(nCP_),
    cpDeltaT0Value_(nCP_),
    cpTimeIndex_(nCP_),
    cpLevel_(nCP_),
    cpAddressing_(nCP_)
{
    init();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ExtendedTime::~ExtendedTime()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
ExtendedTime::loop()
{
    addCP();

    bool running = Time::loop();

    if (debug)
    {
        Info<< "Checkpointing database:" << endl
            << "cp#, addr, lvl, val, idx" << endl;
        forAll(checkpoint_, cpi)
        {
            Info<< cpi << ", "
                << cpAddressing_[cpi] << ", "
                << cpLevel_[cpi] << ", "
                << cpTimeValue_[cpi] << ", "
                << cpTimeIndex_[cpi] << endl;
        }
    }

    if (!running)
    {
        functionObjects().off();

        if (storeCPs_)
        {
            forAll(checkpoint_, cpi)
            {
                forAllConstIter
                (
                    HashTable<regIOobject*>,
                    checkpoint_[cpi],
                    iter
                )
                {
                    iter()->write();
                }
            }
        }
    }

    return running;
}

bool
ExtendedTime::backloop()
{
    bool running = inverseLoop();
    endTime_ = value();

    if (cpIsAvailable())
    {
        restoreCurrentTimeCP();
        deltaT0_ = reverseDeltaTSave_;
        deltaT_ = reverseDeltaT_;
    }
    else if (isValidTimeStepInReverse(value()))
    {
        restoreLatestTimeCP();
    }

    if (debug)
    {
        Info<< "Checkpointing database:" << endl
            << "cp#, addr, lvl, val, idx" << endl;
        forAll(checkpoint_, cpi)
        {
            Info<< cpi << ", "
                << cpAddressing_[cpi] << ", "
                << cpLevel_[cpi] << ", "
                << cpTimeValue_[cpi] << ", "
                << cpTimeIndex_[cpi] << endl;
        }
    }

    if
    (
        !running
     && !writeTime_
    )
    {
        writeNow();
    }

    return running;
}

bool
ExtendedTime::forloop()
{
    if (!cpIsAvailable())
        addCP();

    bool running = Time::loop();

    if (debug)
    {
        Info<< "Checkpointing database:" << endl
            << "cp#, addr, lvl, val, idx" << endl;
        forAll(checkpoint_, cpi)
        {
            Info<< cpi << ", "
                << cpAddressing_[cpi] << ", "
                << cpLevel_[cpi] << ", "
                << cpTimeValue_[cpi] << ", "
                << cpTimeIndex_[cpi] << endl;
        }
    }

    return running;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

ExtendedTime&
ExtendedTime::operator-=(const dimensionedScalar& deltaT)
{
    return operator-=(deltaT.value());
}

ExtendedTime&
ExtendedTime::operator-=(const scalar deltaT)
{
    setDeltaT(deltaT);
    return operator--();
}

ExtendedTime&
ExtendedTime::operator--()
{
    deltaT0_ = reverseDeltaTSave_;
    reverseDeltaTSave_ = reverseDeltaT_;

    // Save old time name
    const word oldTimeName = dimensionedScalar::name();

    setTime(value() - reverseDeltaT_, timeIndex_ - reverseTimestep_);

    if (!subCycling_)
    {
        // If the time is very close to zero reset to zero
        if (mag(value()) < 10*SMALL*primalDeltaT_)
        {
            setTime(0.0, timeIndex_);
        }
    }


    // Check that new time representation differs from old one
    if (dimensionedScalar::name() == oldTimeName)
    {
        int oldPrecision = precision_;
        do
        {
            precision_++;
            setTime(value(), timeIndex());
        }
        while (precision_ < 100 && dimensionedScalar::name() == oldTimeName);

        WarningInFunction
            << "Increased the timePrecision from " << oldPrecision
            << " to " << precision_
            << " to distinguish between timeNames at time " << value()
            << endl;

        if (precision_ == 100 && precision_ != oldPrecision)
        {
            // Reached limit.
            WarningInFunction
                << "Current time name " << dimensionedScalar::name()
                << " is the old as the previous one " << oldTimeName
                << endl
                << "    This might result in overwriting old results."
                << endl;
        }
    }


    if (!subCycling_)
    {
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
                label outputIndex = label
                (
                    ((value() - startTime_) - 0.5*primalDeltaT_)
                  / writeInterval_
                );

                if (outputIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = outputIndex;
                }
            }
            break;

            case wcCpuTime:
            {
                label outputIndex = label
                (
                    returnReduce(elapsedCpuTime(), maxOp<double>())
                  / writeInterval_
                );
                if (outputIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = outputIndex;
                }
            }
            break;

            case wcClockTime:
            {
                label outputIndex = label
                (
                    returnReduce(label(elapsedClockTime()), maxOp<label>())
                  / writeInterval_
                );
                if (outputIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = outputIndex;
                }
            }
            break;

            default:
            break;
        }


        // see if endTime needs adjustment to stop at the next run()/end() check
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

        // Override writeTime if one-shot writing
        if (writeOnce_)
        {
            writeTime_ = true;
            writeOnce_ = false;
        }

    }

    return *this;
}

ExtendedTime&
ExtendedTime::operator--(int)
{
    return operator--();
}

} // End namespace FOAM

// ************************************************************************* //
