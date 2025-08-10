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
    (c) 2016-2019 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "signals/sigFpe.H"
#include "db/error/error.H"
#include "global/JobInfo/JobInfo.H"
#include "include/OSspecific.H"
#include "db/IOstreams/IOstreams.H"

#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
    #include <malloc.h>
#elif defined(sgiN32) || defined(sgiN32Gcc)
    #include <sigfpe.h>
#endif

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigFpe::oldAction_;

bool Foam::sigFpe::sigFpeActive_ = false;

void Foam::sigFpe::fillNan(UList<scalar>& lst)
{
    lst = std::numeric_limits<scalar>::signaling_NaN();
}

bool Foam::sigFpe::mallocNanActive_ = false;


#ifdef LINUX
extern "C"
{
    extern void* __libc_malloc(size_t size);

    // Override the GLIBC malloc to support mallocNan
    void* malloc(size_t size)
    {
        if (Foam::sigFpe::mallocNanActive_)
        {
            return Foam::sigFpe::mallocNan(size);
        }
        else
        {
            return __libc_malloc(size);
        }
    }
}

void* Foam::sigFpe::mallocNan(size_t size)
{
    // Call the low-level GLIBC malloc function
    void * result = __libc_malloc(size);

    // Initialize to signalling NaN
    UList<scalar> lst(reinterpret_cast<scalar*>(result), size/sizeof(scalar));
    sigFpe::fillNan(lst);

    return result;
}
#endif


#if defined(LINUX_GNUC) && defined(RUNTIME_DEBUG_FLAGS)
void Foam::sigFpe::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    set(false);
}


Foam::sigFpe::ignore::ignore()
:
    wasActive_(sigFpeActive_)
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

void Foam::sigFpe::set(const bool verbose)
{
    if (!sigFpeActive_ && (env("FOAM_SIGFPE") || env("FOAM_SIGFPE")))
    {
        bool supported = false;

#ifdef RUNTIME_DEBUG_FLAGS
#ifdef LINUX_GNUC
        supported = true;

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }

        sigFpeActive_ = true;

#elif defined(sgiN32) || defined(sgiN32Gcc)
        supported = true;

        sigfpe_[_DIVZERO].abort=1;
        sigfpe_[_OVERFL].abort=1;
        sigfpe_[_INVALID].abort=1;

        sigfpe_[_DIVZERO].trace=1;
        sigfpe_[_OVERFL].trace=1;
        sigfpe_[_INVALID].trace=1;

        handle_sigfpes
        (
            _ON,
            _EN_DIVZERO
          | _EN_INVALID
          | _EN_OVERFL,
            0,
            _ABORT_ON_ERROR,
            nullptr
        );

        sigFpeActive_ = true;

#endif
#endif // RUNTIME_DEBUG_FLAGS


        if (verbose)
        {
            if (supported)
            {
                Info<< "sigFpe : Enabling floating point exception trapping"
                    << " (FOAM_SIGFPE)." << endl;
            }
            else
            {
#ifdef RUNTIME_DEBUG_FLAGS
                Info<< "sigFpe : Floating point exception trapping"
                    << " - not supported on this platform" << endl;
#else // RUNTIME_DEBUG_FLAGS
                Info<< "sigFpe : Floating point exception trapping"
                    << " - not supported in this build" << endl;
#endif // RUNTIME_DEBUG_FLAGS
            }
        }
    }


    if (env("FOAM_SETNAN") || env("FOAM_SETNAN"))
    {
        #ifdef LINUX
        mallocNanActive_ = true;
        #endif

        if (verbose)
        {
            if (mallocNanActive_)
            {
                Info<< "SetNaN : Initialising allocated memory to NaN"
                    << " (FOAM_SETNAN)." << endl;
            }
            else
            {
                Info<< "SetNaN : Initialise allocated memory to NaN"
                    << " - not supported on this platform" << endl;
            }
        }
    }
}


void Foam::sigFpe::unset(const bool verbose)
{
    #ifdef LINUX_GNUC
    // Reset signal
    if (sigFpeActive_)
    {
        if (verbose)
        {
            Info<< "sigFpe : Disabling floating point exception trapping"
                << endl;
        }

        if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }

        // Reset exception raising
        int oldExcept = fedisableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        if (oldExcept == -1)
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }
        sigFpeActive_ = false;
    }
    #endif

    #ifdef LINUX
    // Disable initialization to NaN
    mallocNanActive_ = false;
    #endif
}


// ************************************************************************* //
