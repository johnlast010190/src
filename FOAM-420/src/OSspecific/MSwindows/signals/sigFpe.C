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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2011 Symscape
    (c) 2016-2018 OpenCFD Ltd.

Notes
    Unlike other parts of OSspecific for MSwindows, this file was taken from ESI
    OpenFOAM v1912.

Class
    sigFpe

\*---------------------------------------------------------------------------*/

#include "signals/sigFpe.H"
#include "db/error/error.H"
#include "global/JobInfo/JobInfo.H"
#include "include/OSspecific.H"
#include "db/IOstreams/IOstreams.H"
#include "primitives/bools/Switch/Switch.H"
#include "containers/Lists/UList/UList.H"

#include <float.h>  // For *fp functions
#include <limits>

// File-local functions
#include "signalMacros.C"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::sigFpe::switchFpe_(Foam::debug::optimisationSwitch("trapFpe", 0));
bool Foam::sigFpe::switchNan_(Foam::debug::optimisationSwitch("setNaN", 0));

bool Foam::sigFpe::sigActive_ = false;
bool Foam::sigFpe::nanActive_ = false;

// Saved old FPE signal trapping setting (file-local variable)
static unsigned int oldFpe_ = 0u;


static void clearFpe()
{
    #ifndef Foam_no_sigFpe
    _clearfp();
    _controlfp(oldFpe_, 0xFFFFFFFF);
    #endif
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Can turn on/off via env variable containing a bool (true|false|on|off ...)
// or by the specified flag
static bool isTrue(const char* envName, bool deflt)
{
    const auto str(Foam::getEnv(envName));

    if (str.size())
    {
        Foam::Switch sw(str, true);  // Silently ignores bad input

        if (sw.valid())
        {
            return sw;
        }
    }

    // Env was not set or did not contain a valid bool value
    return deflt;
    }


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigFpe::sigHandler(int)
{
    resetHandler("SIGFPE", SIGFPE);

    jobInfo.signalEnd();        // Update jobInfo file
    error::printStack(Perr);
    clearFpe();
    ::raise(SIGFPE);            // Throw signal (to old handler)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    set(false);
}


Foam::sigFpe::ignore::ignore()
:
    wasActive_(sigFpe::active())
{
    if (wasActive_)
    {
        sigFpe::unset();
}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    unset(false);
}


Foam::sigFpe::ignore::~ignore()
    {
    restore();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::ignore::restore()
{
    if (wasActive_)
    {
        sigFpe::set();
        }
    wasActive_ = false;
}


bool Foam::sigFpe::requested()
{
    return (isTrue("FOAM_SIGFPE", switchFpe_) || isTrue("FOAM_SIGFPE", switchFpe_));
}


void Foam::sigFpe::set(bool verbose)
{
    if (!sigActive_ && requested())
    {
        #ifdef Foam_no_sigFpe

        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";
            Info<< "- disabled on this platform" << endl;
        }

        #else

        oldFpe_ = _controlfp(0, 0);

        const unsigned int newFpe =
        (
            oldFpe_ & ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW)
        );

        _controlfp(newFpe, _MCW_EM);

        setHandler("SIGFPE", SIGFPE, sigHandler);

        sigActive_ = true;

        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";

            if (sigActive_)
            {
                Info<< "enabled (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "- not supported on this platform" << endl;
            }
        }
        #endif
    }


    nanActive_ = false;
    if (isTrue("FOAM_SETNAN", switchNan_) || isTrue("FOAM_SETNAN", switchNan_))
    {
        if (verbose)
        {
            Info<< "setNaN : Initialise allocated memory to NaN "
                << "- not supported on this platform" << endl;
        }
    }
}


void Foam::sigFpe::unset(bool verbose)
{
    if (sigActive_)
    {
        if (verbose)
        {
            Info<< "sigFpe : Disabling floating point exception trapping"
              << endl;
        }

        sigActive_ = false;

        clearFpe();

        resetHandler("SIGFPE", SIGFPE);
    }

    nanActive_ = false;
}


void Foam::sigFpe::fillNan(UList<scalar>& list)
{
    list = std::numeric_limits<scalar>::signaling_NaN();
}


// ************************************************************************* //
