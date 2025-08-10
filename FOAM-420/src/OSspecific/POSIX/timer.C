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

\*---------------------------------------------------------------------------*/

#include <unistd.h>

#include "db/error/error.H"
#include "timer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(timer, 0);

jmp_buf timer::envAlarm;

struct sigaction timer::oldAction_;

unsigned int timer::oldTimeOut_ = 0;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::timer::signalHandler(int)
{
    if (debug)
    {
        InfoInFunction<< "Timed out. Jumping."
            << endl;
    }
    longjmp(envAlarm, 1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timer::timer(const unsigned int newTimeOut)
:
    newTimeOut_(newTimeOut)
{

    if (newTimeOut > 0)
    {
        // Is singleton since handler is static function
        if (oldTimeOut_ != 0)
        {
            FatalErrorInFunction
                << "timer already used."
                << abort(FatalError);
        }

        // Install alarm signal handler:
        // - do not block any signals while in it
        // - clear list of signals to mask
        struct sigaction newAction;
        newAction.sa_handler = timer::signalHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);

        if (sigaction(SIGALRM, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }

        oldTimeOut_ = ::alarm(newTimeOut);

        if (debug)
        {
            InfoInFunction
                << "Installing timeout " << int(newTimeOut_)
                << " seconds"
                << " (overriding old timeout " << int(oldTimeOut_)
                << ")." << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timer::~timer()
{
    if (newTimeOut_ > 0)
    {
        if (debug)
        {
            InfoInFunction
                << "timeOut=" << int(newTimeOut_)
                << " : resetting timeOut to " << int(oldTimeOut_) << endl;
        }

        // Reset timer
        ::alarm(oldTimeOut_);
        oldTimeOut_ = 0;

        // Restore signal handler
        if (sigaction(SIGALRM, &oldAction_, nullptr) < 0)
        {
            FatalErrorInFunction
                << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
