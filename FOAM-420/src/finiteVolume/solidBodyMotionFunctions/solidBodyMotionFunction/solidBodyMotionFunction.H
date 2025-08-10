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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

Namespace
    Foam::solidBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::solidBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    solidBodyMotionFunction.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionFunction_H
#define solidBodyMotionFunction_H

#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "db/dictionary/dictionary.H"
#include "primitives/septernion/septernion.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionFunction
{
protected:

    // Protected data

        //- Time
        const Time& time_;

        //- Motion dictionary
        dictionary SBMFCoeffs_;

        //- Current timestep
        mutable scalar dt_;

        //- Old timestep
        mutable scalar dtOld_;

        //- Initial phase
        // mutable scalar t0_;

        //- Old time
        mutable scalar tOld_;

        //- Current time index
        mutable label updateIndex_;

        //- Is relative motion? (true to support relative motions: AMR)
        Switch incrementalMotion_;

        //- Coordinate frame name
        word frameName_;

        //- Pointer to coordinate frame
        coordinateFrame* coorFramePtr_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidBodyMotionFunction(const solidBodyMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const solidBodyMotionFunction&);


public:

    //- Runtime type information
    TypeName("solidBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidBodyMotionFunction,
            dictionary,
            (
                const dictionary& SBMFCoeffs,
                const Time& runTime,
                const word& frameName
            ),
            (SBMFCoeffs, runTime, frameName)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidBodyMotionFunction,
            registry,
            (
                const objectRegistry& obr,
                const dictionary& SBMFCoeffs,
                const word& frameName
            ),
            (obr, SBMFCoeffs, frameName)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary, Time and frame name
        solidBodyMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        //- Construct from the objectRegistry, SBMFCoeffs dictionary
        // and frame name
        solidBodyMotionFunction
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const = 0;


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary, Time
        // and frameName
        static autoPtr<solidBodyMotionFunction> New
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        //- Select constructed from the objectRegistry, SBMFCoeffs dictionary
        // and frame name
        static autoPtr<solidBodyMotionFunction> New
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );


    //- Destructor
    virtual ~solidBodyMotionFunction();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const = 0;

        //- Motion velocity at specific time (ds/dt)
        virtual vectorTuple velocity() const = 0;

        //- Motion acceleration at specific time (ds/dt)
        virtual vectorTuple acceleration() const;

        //- Default centre of rotation is coordinate system origin
        virtual const vector& CofR() const;

        //- Default centre of rotation is coordinate system origin
        virtual const vector& CofR0() const;

        //- Force only incremental motion definition
        void setIncrementalMotion(const word& type);

        //- Reset times for relative motions
        void updateTime() const;

        //- Return if the motion is relative or not
        Switch& isIncrementalMotion() const;

        //- Does it allow outer corrector motion?
        virtual bool outerCorrectorMotion() const
        {
            return false;
        }

        //- Does the motion have a coordinate frame?
        bool hasFrame() const
        {
            return coorFramePtr_ != nullptr;
        }

        //- Frame access
        coordinateFrame& frame() const
        {
            if (!coorFramePtr_)
            {
                FatalErrorIn
                (
                    "solidBodyMotionFunction::frame()"
                )
                << "No coordinate frame defined for motion function" << nl
                << abort(FatalError);
            }
            return *coorFramePtr_;
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);

        //- Write in dictionary format
        virtual void writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
