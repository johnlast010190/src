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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrameState/coordinateFrameState.H"
#include "fvMesh/fvMesh.H"
#include "primitives/transform/transform.H"
#include "coordinate/systems/cartesianCS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateFrameState::coordinateFrameState
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& frameName
)
:
    frameStateName_(frameName + "State"),
    isAutonomousDynamic_(false),
    incrementalMotion_(false),
    outerCorrectorMotion_(false),
    frameName_(frameName),
    coordinateSystem_(nullptr),
    coordinateSystemTime0_(nullptr),
    oldUpdateIndex_(obr.time().timeIndex()),
    state0Ptr_(nullptr),
    obr_(obr),
    transformations_(1, septernion::I),
    velocities_(1, vectorTuple(Zero, Zero)),
    accelerations_(1, vectorTuple(Zero, Zero)),
    updateIndex_(-1),
    outerCorrectorUpdate_(1)
{
    read(dict);
}


Foam::coordinateFrameState::coordinateFrameState
(
    const coordinateFrameState& frameState
)
:
    frameStateName_(frameState.frameStateName_),
    isAutonomousDynamic_(frameState.isAutonomousDynamic_),
    incrementalMotion_(frameState.incrementalMotion_),
    outerCorrectorMotion_(frameState.outerCorrectorMotion_),
    frameName_(frameState.frameName_),
    coordinateSystem_(frameState.coordinateSystem_().clone()),
    coordinateSystemTime0_(frameState.coordinateSystemTime0_().clone()),
    oldUpdateIndex_(frameState.oldUpdateIndex_),
    state0Ptr_(nullptr),
    obr_(frameState.obr_),
    transformations_(frameState.transformations_),
    velocities_(frameState.velocities_),
    accelerations_(frameState.accelerations_),
    updateIndex_(frameState.updateIndex_),
    outerCorrectorUpdate_(frameState.outerCorrectorUpdate_)
{
    if (frameState.state0Ptr_)
    {
        state0Ptr_ = new coordinateFrameState(*frameState.state0Ptr_);
    }
}


Foam::coordinateFrameState::coordinateFrameState
(
    const word& newName,
    const coordinateFrameState& frameState
)
:
    frameStateName_(newName),
    isAutonomousDynamic_(frameState.isAutonomousDynamic_),
    incrementalMotion_(frameState.incrementalMotion_),
    outerCorrectorMotion_(frameState.outerCorrectorMotion_),
    frameName_(frameState.frameName_),
    coordinateSystem_(frameState.coordinateSystem_().clone()),
    coordinateSystemTime0_(frameState.coordinateSystemTime0_().clone()),
    oldUpdateIndex_(frameState.oldUpdateIndex_),
    state0Ptr_(nullptr),
    obr_(frameState.obr_),
    transformations_(frameState.transformations_),
    velocities_(frameState.velocities_),
    accelerations_(frameState.accelerations_),
    updateIndex_(frameState.updateIndex_),
    outerCorrectorUpdate_(frameState.outerCorrectorUpdate_)
{
    if (frameState.state0Ptr_)
    {
        state0Ptr_ =
            new coordinateFrameState(newName + "_0", *frameState.state0Ptr_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateFrameState::~coordinateFrameState()
{
    delete state0Ptr_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::septernion& Foam::coordinateFrameState::transformation
(
    bool addParentFrames,
    label nCorr
) const
{
    return transformations_[nCorr != -1 ? nCorr : nFrameCorrector()];
}


void Foam::coordinateFrameState::updateCoordinateSystem() const
{
    // Note the coordinate frame might need to be updated even if it doesn't
    // have it's own motion (isAutonomousDynamic_ == false) because it might be
    // a child of a dynamic frame (e.g. a rotating fan in a moving car)
    vector origin, e1, e2;

    // Compute new origin and axis based on CS at t=0
    if (incrementalMotion_)
    {
        origin = oldTime().coordinateSystem_().origin();
        e1 = oldTime().coordinateSystem_().e1();
        e2 = oldTime().coordinateSystem_().e2();
    }
    else
    {
        origin = coordinateSystemTime0_().origin();
        e1 = coordinateSystemTime0_().e1();
        e2 = coordinateSystemTime0_().e2();
    }

    // Apply transformation (copied from transformPoints)
    // if the transformation is local this isn't needed
    forAll(transformations_, transi)
    {
        const septernion trans = this->transformation(true, transi);
        const vector T = trans.t();
        const tensor R = trans.r().R();

        origin = transform(R, origin - T);

        if (mag(R - I) > SMALL)
        {
            // Translation + Rotation
            // Vectors are only rotated!
            e1 = transform(R, e1);
            e2 = transform(R, e2);
        }
    }

    coordinateSystem_().update(origin, e1, e2);
}


const Foam::coordinateSystem& Foam::coordinateFrameState::coorSys() const
{
    return coordinateSystem_();
}


const Foam::coordinateSystem& Foam::coordinateFrameState::coorSys0() const
{
    return coordinateSystemTime0_();
}


const Foam::septernion&
Foam::coordinateFrameState::decoupledTransformation(label nCorr) const
{
    return transformations_[nCorr != -1 ? nCorr : nFrameCorrector()];
}


Foam::vector Foam::coordinateFrameState::axis() const
{
    return coorSys().e3();
}


Foam::vector Foam::coordinateFrameState::axis0() const
{
    return coorSys0().e3();
}


const Foam::vectorTuple& Foam::coordinateFrameState::velocity() const
{
    return velocities_[nFrameCorrector()];
}


const Foam::vectorTuple& Foam::coordinateFrameState::acceleration() const
{
    return accelerations_[nFrameCorrector()];
}


const Foam::word& Foam::coordinateFrameState::frameStateName() const
{
    return frameStateName_;
}


bool Foam::coordinateFrameState::isDynamic() const
{
    return isAutonomousDynamic_;
}


void Foam::coordinateFrameState::resetDynamic(bool isDynamic)
{
    isAutonomousDynamic_ = isDynamic;
}


void Foam::coordinateFrameState::resetUpdate(label outerCorrectorUpdate) const
{
    updateIndex_ = -1;
    outerCorrectorUpdate_ = outerCorrectorUpdate;
}


Foam::label Foam::coordinateFrameState::nFrameCorrector() const
{
    return outerCorrectorUpdate_ - 1;
}


bool& Foam::coordinateFrameState::outerCorrectorMotion() const
{
    return outerCorrectorMotion_;
}


Foam::Switch& Foam::coordinateFrameState::isIncrementalMotion() const
{
    return incrementalMotion_;
}


bool Foam::coordinateFrameState::isUpdated() const
{
    if (updateIndex_ != obr_.time().timeIndex())
    {
        return false;
    }

    return true;
}


const Foam::coordinateFrameState& Foam::coordinateFrameState::oldTime() const
{
    if (!state0Ptr_)
    {
        state0Ptr_ =
            new coordinateFrameState(this->frameStateName_ + "_0", *this);
    }
    else
    {
        storeOldTimes();
    }

    return *state0Ptr_;
}


void Foam::coordinateFrameState::storeOldTimes() const
{
    if
    (
        state0Ptr_
     && oldUpdateIndex_ != obr_.time().timeIndex()
     &&
       !(
            frameStateName_.size() > 2
         && frameStateName_.endsWith("_0")
        )
    )
    {
        storeOldTime();
    }

    // Correct time index
    oldUpdateIndex_ = obr_.time().timeIndex();
}


void Foam::coordinateFrameState::storeOldTime() const
{
    if (state0Ptr_)
    {
        state0Ptr_->storeOldTime();

        state0Ptr_->forceAssign(*this);
        (*state0Ptr_).oldUpdateIndex_ = oldUpdateIndex_;
    }
}


Foam::label Foam::coordinateFrameState::nOldTimes() const
{
    if (state0Ptr_)
    {
        return state0Ptr_->nOldTimes() + 1;
    }
    else
    {
        return 0;
    }
}


void Foam::coordinateFrameState::readStateData
(
    coordinateFrameState& state,
    const dictionary& dict
)
{
    // Backward compatibility
    if (dict.found("transformation"))
    {
        state.transformations_ =
            List<septernion>(dict.lookup<septernion>("transformation"));
        state.velocities_ =
            List<vectorTuple>(dict.lookup<vectorTuple>("velocity"));
        state.accelerations_ =
            List<vectorTuple>(dict.lookup<vectorTuple>("acceleration"));
    }
    else
    {
        state.transformations_ =
            dict.lookup<List<septernion>>("transformations");
        state.velocities_ = dict.lookup<List<vectorTuple>>("velocities");
        state.accelerations_ = dict.lookup<List<vectorTuple>>("accelerations");
    }
}


bool Foam::coordinateFrameState::read(const dictionary& dict)
{
    coordinateSystemTime0_.reset(frameData::frameToGlobal(dict, obr_));

    IOdictionary stateDictAll
    (
        IOobject
        (
            "referenceFramesState",
            obr_.parent().time().timeName(),
            "uniform",
            obr_.parent(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if
    (
        stateDictAll.typeHeaderOk<dictionary>(true)
     && stateDictAll.found(frameName_)
    )
    {
        // Making copy of the dict to get rid of binary/ASCII issues
        const dictionary stateDict(stateDictAll.subDict(frameName_));

        readStateData(*this, stateDict);

        coordinateSystem_.reset
        (
            coordinateSystem::New(stateDict, coordinateSystem::typeName)
        );

        if (stateDict.found("oldCoordinateFrame"))
        {
            coordinateFrameState& oldState =
                const_cast<coordinateFrameState&>(oldTime());

            // Making copy of the dict to get rid of binary/ASCII issues
            const dictionary oldDict(stateDict.subDict("oldCoordinateFrame"));

            readStateData(oldState, oldDict);

            oldState.coordinateSystem_.reset
            (
                coordinateSystem::New(oldDict, coordinateSystem::typeName)
            );
        }
    }
    else
    {
        // State information wasn't found initialising with time0 coor. system
        coordinateSystem_.reset(coordinateSystemTime0_->clone());
    }

    return true;
}


void Foam::coordinateFrameState::addFrameState
(
    dictionary& dict,
    const coordinateFrameState& state,
    const word& dictName
) const
{
    dictionary frame;
    frame.add("transformations", state.transformations_);
    frame.add("velocities", state.velocities_);
    frame.add("accelerations", state.accelerations_);

    // Write data based on type of coordinate system
    dictionary coorSys;
    coorSys.add("type", state.coorSys().type());
    coorSys.add("e1", state.coorSys().e1());
    coorSys.add("e2", state.coorSys().e2());
    coorSys.add("origin", state.coorSys().origin());
    if (state.coorSys().coordinateTransforms().isLocal())
    {
        coorSys.add("definedInFrame", "parent");
    }
    if (state.coorSys().coordinateTransforms().isLocalOnRead())
    {
        coorSys.add("originallyLocal", "true");
    }
    frame.set("coordinateSystem", coorSys);

    // Add the dict to the main dict
    dict.set(dictName, frame);
}


void Foam::coordinateFrameState::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "referenceFramesState",
            obr_.parent().time().timeName(),
            "uniform",
            obr_.parent(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (dict.found(frameName()))
    {
        dict.remove(frameName());
    }

    addFrameState(dict, *this, frameName());

    // So far one old state is written down
    if (nOldTimes() != 0)
    {
        addFrameState
        (
            dict.subDict(frameName()),
            oldTime(),
            "oldCoordinateFrame"
        );
    }

    // Write the dictionary in ASCII with no compression
    dict.regIOobject::writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::compressionType::UNCOMPRESSED,
        true
    );
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::coordinateFrameState::forceAssign
(
    const coordinateFrameState& cfs
)
{
    // Note for the old time consistency on update the state name isn't updated
    isAutonomousDynamic_ = cfs.isAutonomousDynamic_;
    incrementalMotion_ = cfs.incrementalMotion_;
    frameName_ = cfs.frameName_;
    transformations_ = cfs.transformations_;

    coordinateSystem_.reset(cfs.coordinateSystem_().clone());
    coordinateSystemTime0_.reset(cfs.coordinateSystemTime0_().clone());

    velocities_ = cfs.velocities_;
    accelerations_ = cfs.accelerations_;
    oldUpdateIndex_ = cfs.oldUpdateIndex_;
    updateIndex_= cfs.updateIndex_;
    outerCorrectorMotion_ = cfs.outerCorrectorMotion_;
    outerCorrectorUpdate_ = cfs.outerCorrectorUpdate_;
}


// ************************************************************************* //
