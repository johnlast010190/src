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
    Foam::sigWriteNow

Description
    Signal handler for interupt defined by OptimisationSwitches::writeNowSignal

    Write once and continue.

SourceFiles
    sigWriteNow.C

\*---------------------------------------------------------------------------*/

#ifndef sigWriteNow_H
#define sigWriteNow_H

#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                           Class sigWriteNow Declaration
\*---------------------------------------------------------------------------*/

class sigWriteNow
{
    // Private data

        //- Number of signal to use
        static int signal_;

        //- Pointer to Time
        static Time* runTimePtr_;

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;


    // Private Member Functions

        static void sigHandler(int);


public:

    friend class sigStopAtWriteNow;

    //- wip. Have setter have access to signal_
    friend class addwriteNowSignalToOpt;

    // Constructors

        //- Construct null
        sigWriteNow();

        //- Construct from components
        sigWriteNow(const bool verbose, Time& runTime);


    //- Destructor
    ~sigWriteNow();


    // Member functions

        //- (re)set signal catcher
        static void set(const bool verbose);

        //- Is active?
        bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
