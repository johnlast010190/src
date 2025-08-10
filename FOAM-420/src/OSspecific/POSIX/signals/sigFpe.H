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
    (c) 2018-2019 OpenCFD Ltd.

Class
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Controlled by two env vars:
      - \par FOAM_SIGFPE or FOAM_SIGFPE
        Exception trapping
      - \par FOAM_SETNAN or FOAM_SETNAN
        Initialization of all malloced memory to NaN. If FOAM_SIGFPE or
        FOAM_SIGFPE also set, this will cause usage of uninitialized scalars
        to trigger an abort.

    Can be used either directly through the static member functions or
    through the scope of the object (constructor sets trapping; destructor
    restores original).

SourceFiles
    sigFpe.C

\*---------------------------------------------------------------------------*/

#ifndef sigFpe_H
#define sigFpe_H

#include <signal.h>

#if defined(linux) || defined(linux64) || defined(linuxIA64) || \
    defined(linuxARM7) || defined(linuxPPC64) || defined(linuxPPC64le)
    #define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
    #define LINUX_GNUC
#endif

#include "containers/Lists/UList/UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private data

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

        //- Flag to indicate floating point trapping is enabled
        static bool sigFpeActive_;


    // Static data members

        #if defined(LINUX_GNUC) && defined (RUNTIME_DEBUG_FLAGS)
        //- Handler for caught signals
        static void sigHandler(int);
        #endif


public:


    // Constructors

        //- Construct null
        sigFpe();


    //- Destructor
    ~sigFpe();


    // Member functions

        //- Activate SIGFPE signal handler when FOAM_SIGFPE or FOAM_SIGFPE is %set
        //  Fill memory with NaN when FOAM_SETNAN or FOAM_SETNAN is %set
        static void set(const bool verbose = false);

        //- Deactivate SIGFPE signal handler and NaN memory initialisation
        static void unset(const bool verbose = false);

        //- Flag to indicate mallocNan is enabled
        static bool mallocNanActive_;

        #ifdef LINUX
        //- Malloc function which initializes to NaN
        static void* mallocNan(size_t size);
        #endif

        //- Fill block of data with NaN
        static void fillNan(UList<scalar>&);



    // Helper classes

        //- Helper to locally ignore SIGFPE handling.
        //  Restores the original state of the SIGFPE handler on destruction.
        class ignore
        {
            //- The signal handler state when entering
            bool wasActive_;

            //- No copy construct
            ignore(const ignore&) = delete;

            //- No copy assignment
            void operator=(const ignore&) = delete;

            //- No move construct
            ignore(ignore&&) = delete;

            //- No move assignment
            void operator=(ignore&&) = delete;


        public:

            //- Constructor deactivates any previously active SIGFPE handler
            ignore();

            //- Destructor restores the original state of SIGFPE handler
            ~ignore();

            //- Restore the original state of SIGFPE handler
            void restore();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
