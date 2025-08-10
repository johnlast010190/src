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

Class
    Foam::coordinateFrameState

Description
    Contains state of frame of reference. Stored information are
    transformations and local frame velocities (translation + angular)

\*---------------------------------------------------------------------------*/

#ifndef coordinateFrameState_H
#define coordinateFrameState_H

#include "primitives/septernion/septernion.H"
#include "coordinate/systems/coordinateSystem.H"
#include "db/objectRegistry/objectRegistry.H"
#include "db/Time/Time.H"
#include "primitives/Tuple2/Tuple2.H"
#include "include/demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coordinateFrameState;

/*---------------------------------------------------------------------------*\
                   Class coordinateFrameState Declaration
\*---------------------------------------------------------------------------*/

class coordinateFrameState
{
    // Private Data

        //- Name of the frame state
        word frameStateName_;

        //- If the flag is true the frame can move with its own motion
        // and parent motions. If the flag is false than it can move only
        // with parent motions.
        bool isAutonomousDynamic_;

        //- Motion description based on time0
        mutable Switch incrementalMotion_;

        //- Does motion description support outer corrector increments?
        mutable bool outerCorrectorMotion_;

        //- CoordinateFrame name
        word frameName_;

        //- Current coordinate system
        mutable autoPtr<coordinateSystem> coordinateSystem_;

        //- Coordinate system at time 0
        mutable autoPtr<coordinateSystem> coordinateSystemTime0_;

        //- Update index for update of the old-time value
        mutable label oldUpdateIndex_;

        //- Pointer to old time field
        mutable coordinateFrameState* state0Ptr_;


protected:

    // Protected Data Members

        //- Registry access
        const objectRegistry& obr_;

        //- Local motion function transformation (displacement information)
        mutable List<septernion> transformations_;

        //- Current velocity (translational, rotational)
        mutable List<vectorTuple> velocities_;

        //- Current acceleration (translational, rotational)
        mutable List<vectorTuple> accelerations_;

        //- Update index
        mutable label updateIndex_;

        //- Outer corrector frame update
        mutable label outerCorrectorUpdate_;


public:


    //- Constructors
        coordinateFrameState
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& frameName
        );

        //- Construct as copy
        coordinateFrameState(const coordinateFrameState& frameState);

        //- Construct as copy resetting name
        coordinateFrameState
        (
            const word& newName,
            const coordinateFrameState& frameState
        );

    //- Destructor
        virtual ~coordinateFrameState();


    // Member Functions

        //- Update coordinate system with the transformation septernion
        void updateCoordinateSystem() const;

        //- Store the old-time fields
        void storeOldTimes() const;

        //- Store the old-time field
        void storeOldTime() const;

        //- Return the number of old time states stored
        label nOldTimes() const;

        //- Return old time state
        const coordinateFrameState& oldTime() const;


    // Member operators

        void forceAssign(const coordinateFrameState&);


    // Access Functions

        //- Return name of current frame
        const word& frameName() const
        {
            return frameName_;
        }

        //- Coordinate systems (if system is static return always time 0)
        const coordinateSystem& coorSys() const;

        //- Reference to old coordinate system
        const coordinateSystem& coorSys0() const;

        //- Return the frame motion transformation septernion
        const septernion& decoupledTransformation(label nCorr = -1) const;

        //- Return the frame motion transformation septernion
        virtual const septernion& transformation
        (
            bool addParentFrames = true,
            label nCorr = -1
        ) const;

        //- Axis
        vector axis() const;

        //- Initial axis
        vector axis0() const;

        //- Translational, rotational velocity (should be septernion?)
        const vectorTuple& velocity() const;

        //- Translational, rotational acceleration (should be septernion?)
        const vectorTuple& acceleration() const;

        //- Object name
        const word& frameStateName() const;

        //- Return if the frame is dynamic
        bool isDynamic() const;

        //- Set if the frame is dynamic
        void resetDynamic(bool isDynamic);

        //- Reset this time update
        void resetUpdate(label outerCorrectorUpdate = 1) const;

        //- Returns number that corresponds to outer corrector number when it
        //  switched on and the motion supports outer corrector sequences.
        //  Otherwise it returns 1.
        label nFrameCorrector() const;

        //- Does the motion supports outer corrector sequences?
        bool& outerCorrectorMotion() const;

        //- Return if the it is time zero motion
        Switch& isIncrementalMotion() const;

        //- Access to update systems
        bool isUpdated() const;

        //- Backward compatible read properties from given dictionary
        void readStateData
        (
            coordinateFrameState& state,
            const dictionary& dict
        );

        //- Update properties from given dictionary
        bool read(const dictionary& dict);

        //- Support function to create dictionary to write
        void addFrameState
        (
            dictionary& dict,
            const coordinateFrameState& state,
            const word& dictName
        ) const;

        //- Write to stream
        void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//end namespace FOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif


// ************************************************************************* //
