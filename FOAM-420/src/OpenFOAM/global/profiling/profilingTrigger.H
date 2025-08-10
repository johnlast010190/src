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
    (c) 2009-2016 Bernhard Gschaider
    (c) 2016-2107 OpenCFD Ltd.

Class
    Foam::profilingTrigger

Description
    Triggers for starting/stopping code profiling.

SourceFiles
    profilingTrigger.C

\*---------------------------------------------------------------------------*/

#ifndef profilingTrigger_H
#define profilingTrigger_H

#include "clockTime/clockTime.H"
#include "primitives/strings/string/string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class profilingInformation;

/*---------------------------------------------------------------------------*\
                      Class profilingTrigger Declaration
\*---------------------------------------------------------------------------*/

class profilingTrigger
{
    // Private Data Members

        //- The timer for the profiling information
        clockTime clock_;

        //- The profiling information
        profilingInformation *ptr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        profilingTrigger(const profilingTrigger&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const profilingTrigger&) = delete;


public:

    // Constructors

        //- Construct profiling with given description.
        //  Descriptions beginning with 'application::' are reserved for
        //  internal use.
        profilingTrigger(const char* name);

        //- Construct profiling with given description.
        //  Descriptions beginning with 'application::' are reserved for
        //  internal use.
        profilingTrigger(const string& name);


    //- Destructor
    ~profilingTrigger();


    // Member Functions

      // Access

        //- True if the triggered profiling is active
        bool running() const;


      // Edit

        //- Stop triggered profiling
        void stop();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macros

//- Define profiling with specified name and description string
//  This is required if the description contains space, colons etc.
//  \sa addProfiling0
//  \sa endProfiling
#define addProfiling(name,descr)                                               \
    ::Foam::profilingTrigger  profilingTriggerFor##name(descr)

//- Define profiling with specified name and description correspond to the name
//  \sa addProfiling
//  \sa endProfiling
#define addProfiling0(name)                                                    \
    ::Foam::profilingTrigger  profilingTriggerFor##name(#name)

//- Define profiling with specified name and description correspond to the
//  compiler-defined function name string:
//  \sa addProfiling
//  \sa endProfiling
#ifdef __GNUC__
    #define addProfilingInFunction(name)                                       \
    ::Foam::profilingTrigger  profilingTriggerFor##name(__PRETTY_FUNCTION__)
#else
    #define addProfilingInFunction(name)                                       \
    ::Foam::profilingTrigger  profilingTriggerFor##name(__func__)
#endif

//- Remove profiling with specified name
//  \sa addProfiling
//  \sa addProfiling0
#define endProfiling(name)    profilingTriggerFor##name.stop()


#endif

// ************************************************************************* //
