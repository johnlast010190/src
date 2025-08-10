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
    (c) 2019 OpenCFD Ltd.

Description
    File-local code for setting/resetting signal handlers.

SourceFiles
    signalMacros.C

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include <csignal>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Saved old signal trapping setting (file-local variable)
static __p_sig_fn_t oldAction_ = SIG_DFL;


static void resetHandler(const char *what, int sigNum)
{
    const __p_sig_fn_t prev = ::signal(sigNum, oldAction_);
    oldAction_ = SIG_DFL;

    if (SIG_ERR == prev)
    {
        FatalError
            << "Cannot unset " << what << " signal (" << sigNum
            << ") trapping" << endl
            << abort(FatalError);
    }
}


static void setHandler(const char *what, int sigNum, void (*handler)(int))
{
    oldAction_ = ::signal(sigNum, handler);

    if (SIG_ERR == oldAction_)
    {
        FatalError
            << "Could not set " << what << " signal (" << sigNum
            << ") trapping" << endl
            << abort(FatalError);
    }
}

} // End namespace Foam


// ************************************************************************* //
