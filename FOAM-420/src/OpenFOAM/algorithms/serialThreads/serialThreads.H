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
    (c) 2019 Esi Ltd.

Description
    Run threads in sequence, pausing execution and switching between them at
    an arbitrary point in the code. Allows load-balancing between regions.
     - create() creates a function with threads
     - run() runs the function and returns when all threads have completed
     - yeildToNext() pauses the currently running thread and switches to the
    next one.
     - level() waits until all threads reach this point in the code
     - pauseSwitching/resumeSwitching() suspends or resumes thread switching
    when yieldToNext is called.  Typically used to surround a section of code
    that is not load-balanced between regions - e.g. fvOptions which are not
    the same for different regions and so may contain unmatched parallel
    synchronisation calls.

\*---------------------------------------------------------------------------*/

#ifndef serialThreads_H
#define serialThreads_H

#include "primitives/ints/label/label.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "algorithms/DelayedFunctionCall/DelayedFunctionCall.H"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class serialThreads
{
    protected:

        static std::vector<std::thread> threads_;

        static label activeThread_;
        static std::mutex activeThreadMutex_;
        static std::condition_variable masterNotify_;
        static std::vector<std::condition_variable> threadNotify_;

        // Total number of threads in the pool
        static label numThreads_;
        // Number of threads actually activated (running)
        static label numActiveThreads_;
        static thread_local label thisThreadNum_;

        // Records the number of threads that have paused; only resumes when
        // the same number have un-paused
        static label paused_;

        // Flag set to tell threads to exit
        static bool quit_;

        static autoPtr<DelayedFunctionCall<void>> functionCall_;

    protected:
        static void wait();
        static void activateNext();
        static void activateMaster();
        static void activateThread(const label iThread);

        static void threadEntryPoint(const label nThread);

    public:
        // Sets up thread pool (can involve significant overhead)
        static void create(const label maxNumThreads);

        // Runs the member function of ObjType in threads, sequentially.
        // numThreads must be less than or equal to maxNumThreads in create()
        template<class ObjType, typename... ArgsType>
        static void run
        (
            const label numThreads,
            ObjType& obj,
            void (ObjType::*memFuncPtr)(ArgsType...),
            ArgsType&&... args
        );

        // Returns whether currently active (run() called)
        static bool active();

        // Ends the threads
        static void end();

        // Pause the currently running thread and switch to the next one
        // in sequence. Typically put before a parallel synchronisation point
        // where work across the different threads needs to load-balance.
        static void yieldToNext();

        // Pause currently running thread and switch to master
        static void yieldToMaster();

        // Return the number (from 0 to numThreads) of the present thread
        static label thisThreadNum();

        // Return number of threads currently running
        static label numThreads();

        // Hold threads as necessary to get them back in the original order -
        // typically used after code which is not the same between different
        // regions, but allows for partial load-balancing across regions.
        static void level();

        // Disable further calls to yieldToNext, typically put before a section
        // of code that is not the same between different regions. As a result,
        // no load balancing occurs across regions (rather use 'level()').
        static void pauseSwitching();

        // Re-enable calls to yeildToNext
        static void resumeSwitching();
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
