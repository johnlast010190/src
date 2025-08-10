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

Class
    Foam::timeControl

Description
    General time dependent execution controller.
    The default to execute every time-step.

SourceFiles
    timeControl.C

\*---------------------------------------------------------------------------*/

#ifndef timeControl_H
#define timeControl_H

#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeControl Declaration
\*---------------------------------------------------------------------------*/

class timeControl
{
public:

    //- The time control options
    enum timeControls
    {
        ocTimeStep,           //!< execution is coupled to the time-step
        ocWriteTime,          //!< execution is coupled to the write-time
        ocOutputTime,         //!< execution is coupled to the output-time
        ocOutputTimeAndEnd,   //!< execution at output-time and end of run
        ocAdjustableRunTime,  //!< Adjust time step for execution
        ocRunTime,            //!< run time for execution
        ocClockTime,          //!< clock time for execution
        ocCpuTime,            //!< CPU time for execution
        ocOnEnd,              //!< on end of run
        ocNone                //!< no execution
    };


private:

    // Private data

        //- Time object
        const Time& time_;

        //- Prefix
        const word prefix_;

        //- String representation of timeControls enums
        static const NamedEnum<timeControls, 10> timeControlNames_;

        //- Type of time control
        timeControls timeControl_;

        //- Execution interval steps for timeStep mode
        //  a value <= 1 means execute at every time step
        label intervalSteps_;

        //- Execution interval
        scalar interval_;

        //- Index of previous execution
        label executionIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        timeControl(const timeControl&);

        //- Disallow default bitwise assignment
        void operator=(const timeControl&);


public:

    // Constructors

        //- Construct from Time object and dictionary
        timeControl
        (
            const Time&,
            const dictionary&,
            const word& prefix
        );


    //- Destructor
    ~timeControl();


    // Member Functions

        //- Helper function to identify if a timeControl object is present
        //  in the dictionary
        static bool entriesPresent(const dictionary& dict, const word& prefix);

        //- Read from dictionary
        void read(const dictionary&);

        //- Return Time
        inline const Time& time() const;

        //- Flag to indicate whether to execute
        bool execute();

        //- Return control
        inline timeControls control() const;

        //- Return interval
        inline scalar interval() const;

        //- Return the index of the previous execution
        inline label executionIndex() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "db/functionObjects/timeControl/timeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
