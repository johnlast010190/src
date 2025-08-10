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
    (c) 2016 blueCAPE Ltda.


Modifications
    This file is based on the original version for POSIX:
        OpenFOAM/src/OSspecific/POSIX/

    This file has been created by blueCAPE's unofficial mingw patches for
    OpenFOAM.
    For more information about these patches, visit:
        http://bluecfd.com/Core

    Modifications made:
      - Derived from the patches for blueCFD 2.1 and 2.2.
      - Adapted for blueCFD-Core 2017.
      - Minor changes for FOAM CORE

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "include/OSspecific.H"

//Undefine this macro, otherwise it will collide with Windows' definitions
#undef DebugInformation

#include "stack_trace.h"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void error::safePrintStack(std::ostream& os)
{
    std::stringstream callstacktext
    (
        std::stringstream::in | std::stringstream::out
    );

    os << "Generating stack trace..." << std::endl;

    StackTrace *traceUs = new StackTrace();

    if (traceUs!=nullptr)
    {
        traceUs->OutputToStream(&callstacktext);

        delete traceUs;
        traceUs=nullptr;

        os << callstacktext.str().data();
    }
    else
    {
        os << "We're sorry, but the application crashed and safe stack tracing "
              "isn't workign in this current windows application in FOAM "
              "please contact ESI support for more information."
           << endl;
    }
}

void error::printStack(Ostream& os)
{
    std::stringstream callstacktext
    (
        std::stringstream::in | std::stringstream::out
    );

    os << "Generating stack trace..." << endl;

    StackTrace *traceUs = new StackTrace();

    if (traceUs!=nullptr)
    {
        traceUs->OutputToStream(&callstacktext);

        delete traceUs;
        traceUs=nullptr;

        os << callstacktext.str().data();
    }
    else
    {
        os << "We're sorry, but the application crashed and safe stack tracing "
              "isn't workign in this current windows application in FOAM "
              "please contact ESI support for more information."
           << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
