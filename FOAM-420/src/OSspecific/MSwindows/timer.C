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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2011 Symscape

Modifications
    This file is based on the original version for POSIX:
        OpenFOAM/src/OSspecific/POSIX/

    This file was developed for Windows by:
        Copyright            : (C) 2011 Symscape
        Website              : www.symscape.com

    This copy of this file has been created by blueCAPE's unofficial mingw
    patches for OpenFOAM.
    For more information about these patches, visit:
        http://bluecfd.com/Core

    Modifications made:
      - Derived from the patches for blueCFD 2.1 and 2.2.
      - Adjusted the code to OpenFOAM 2.2.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "MSwindows.H"
#include "timer.H"

// Undefine DebugInformation, because we don't use it here and it collides with a
// macro in windows.h
#undef DebugInformation

#ifndef WINVER
#define WINVER 0x0500 // To access CreateTimerQueueTimer
#else
#if (WINVER < 0x0500)
#undef WINVER
#define WINVER 0x0500 // To access CreateTimerQueueTimer
#endif
#endif

#include <windows.h>

#define SIGALRM 14


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(timer, 0);

jmp_buf timer::envAlarm;

__p_sig_fn_t timer::oldAction_ = SIG_DFL;

static HANDLE hTimer_ = nullptr;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::timer::signalHandler(int)
{
    if (debug)
    {
        InfoInFunction
            << "timed out. Jumping."
            << endl;
    }
    ::longjmp(envAlarm, 1);
}


static VOID CALLBACK timerExpired(PVOID lpParam, BOOLEAN TimerOrWaitFired)
{
    ::raise(SIGALRM);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
Foam::timer::timer(const unsigned int newTimeOut)
:
    newTimeOut_(newTimeOut)
{

    if (newTimeOut > 0)
    {
        // Is singleton since handler is static function
        if (nullptr != hTimer_)
        {
            FatalErrorInFunction
                << "timer already used."
                << abort(FatalError);
        }

        // Install alarm signal handler:
        oldAction_ = ::signal(SIGALRM, &Foam::timer::signalHandler);

        if (SIG_ERR == oldAction_)
        {
            oldAction_ = SIG_DFL;

            FatalErrorInFunction
                << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }

        if (debug)
        {
            Info<< "Foam::timer::timer(const unsigned int) : "
                << " installing timeout " << int(newTimeOut_)
                << " seconds." << endl;
        }

        const bool success =
          ::CreateTimerQueueTimer(&hTimer_,
                                  nullptr,
                                  WAITORTIMERCALLBACK(timerExpired),
                                  nullptr ,
                                  newTimeOut * 1000,
                                  0, 0);

        if (!success)
        {
            hTimer_ = nullptr;
            FatalErrorInFunction
                << "CreateTimerQueueTimer, "
                << MSwindows::getLastError()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timer::~timer()
{
    if (newTimeOut_ > 0)
    {
        // Reset timer
        const bool timerSuccess =
          ::DeleteTimerQueueTimer(nullptr, hTimer_, nullptr);
        hTimer_ = nullptr;

        if (!timerSuccess)
        {
            FatalErrorInFunction
                << "DeleteTimerQueueTimer, "
                << MSwindows::getLastError()
                << abort(FatalError);
        }

        if (debug)
        {
            InfoInFunction
                << "timeOut=" << int(newTimeOut_) << endl;
        }

        const __p_sig_fn_t signalSuccess = signal(SIGALRM, oldAction_);
        oldAction_ = SIG_DFL;

        // Restore signal handler
        if (SIG_ERR == signalSuccess)
        {
            FatalErrorInFunction
                << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
