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
    Foam::sigSegv

Description
    Signal handler for SEGV interupt.

    The standard interupt handler is overridden to ensure that the
    runningJob file is removed.
    Can be used either directly through the static member functions or
    through the scope of the object (constructor sets trapping; destructor
    restores original).

See also
    Foam::JobInfo

SourceFiles
    sigSegv.C

\*---------------------------------------------------------------------------*/

#ifndef sigSegv_H
#define sigSegv_H

#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigSegv Declaration
\*---------------------------------------------------------------------------*/

class sigSegv
{
    // Private data

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

        //- Flag to indicate signal trapping is enabled
        static bool sigActive_;


    // Private Member Functions

        //- Handler for caught signals
        static void sigHandler(int);


public:

    // Constructors

        //- Construct null
        sigSegv();


    //- Destructor
    ~sigSegv();


    // Member functions

        //- Activate SIGSEGV signal handler
        static void set(const bool verbose);

        //- Deactivate SIGSEGV signal handler
        static void unset(const bool verbose);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
