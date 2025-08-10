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
    (c) 2019-2022 Esi Ltd.

Application
    serialThreads.C

\*---------------------------------------------------------------------------*/

#include "algorithms/serialThreads/serialThreads.H"
#include "db/error/error.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * Static member initialisation  * * * * * * * * * * * * * //

std::vector<std::thread> Foam::serialThreads::threads_;
Foam::label Foam::serialThreads::activeThread_ = -1;
std::mutex Foam::serialThreads::activeThreadMutex_;
std::condition_variable Foam::serialThreads::masterNotify_;
std::vector<std::condition_variable> Foam::serialThreads::threadNotify_;
Foam::label Foam::serialThreads::numThreads_ = 0;
Foam::label Foam::serialThreads::numActiveThreads_ = 0;
thread_local Foam::label Foam::serialThreads::thisThreadNum_ = -1;
Foam::label Foam::serialThreads::paused_ = 0;
bool Foam::serialThreads::quit_ = false;
Foam::autoPtr<Foam::DelayedFunctionCall<void>> Foam::serialThreads::functionCall_;

// * * * * * * * * * Protected static functions  * * * * * * * * * * * * * * //

void Foam::serialThreads::wait()
{
    // Put this thread into a wait; checks whether active before releasing
    {
        std::unique_lock<std::mutex> waitLock(activeThreadMutex_);
        // Note that the mutex is released while
        // thread is waiting, then reacquired
        if (thisThreadNum_ < 0)
        {
            masterNotify_.wait
            (
                waitLock,
                []{return thisThreadNum_ == activeThread_;}
            );
        }
        else
        {
            threadNotify_[thisThreadNum_].wait
            (
                waitLock,
                []{return thisThreadNum_ == activeThread_;}
            );
        }
    }
}

void Foam::serialThreads::activateNext()
{
    // Activate the next thread
    {
        std::lock_guard<std::mutex> modifyLock(activeThreadMutex_);
        activeThread_ = (activeThread_+1) % numActiveThreads_;
        threadNotify_[activeThread_].notify_one();
    }
}

void Foam::serialThreads::activateMaster()
{
    // Activate the master thread
    {
        std::lock_guard<std::mutex> modifyLock(activeThreadMutex_);
        activeThread_ = -1;
        masterNotify_.notify_one();
    }
}

void Foam::serialThreads::activateThread(const label iThread)
{
    // Activate the numbered thread
    {
        std::lock_guard<std::mutex> modifyLock(activeThreadMutex_);
        activeThread_ = iThread;
        if (iThread < 0)
        {
            masterNotify_.notify_one();
        }
        else
        {
            threadNotify_[activeThread_].notify_one();
        }
    }
}


void Foam::serialThreads::threadEntryPoint(const label nThread)
{
    // Entry point for threads

    thisThreadNum_ = nThread;

    // Indicate that we are ready to create the next thread
    yieldToMaster();

    while (!quit_)
    {
        if (functionCall_.valid())
        {
            functionCall_->call();
        }
        else
        {
            FatalErrorInFunction
                << "Internal error: Serial thread activated without "
                << "a function call." << endl
                << abort(FatalError);
        }

        if (paused_)
        {
            FatalErrorInFunction
                << "Mismatch in number of pause/resume calls: "
                << "Threads were not resumed before ending."
                << abort(FatalError);
        }

        // Check that all threads are in sync and wait until they are

        static List<bool> threadsArrived(numThreads_, false);
        static bool releasing = false;
        threadsArrived[thisThreadNum_] = true;
        bool lastArrival = true;
        for (label threadi = 0; threadi < numActiveThreads_; threadi++)
        {
            lastArrival &= threadsArrived[threadi];
            bool outOfOrder =
                (!threadsArrived[threadi] && threadi != thisThreadNum_+1);
            bool waited(false);
            if (!threadsArrived[threadi])
            {
                activateThread(threadi);
                wait();
                waited = true;
            }
            if (outOfOrder || (!threadsArrived[threadi] && !releasing))
            {
                WarningInFunction << nl
                    << "Threads did not complete in order." << nl
                    << "This probably means that different regions in a " << nl
                    << "consolidated group are not using the same schemes."
                    << nl << endl;
            }
            while (!threadsArrived[threadi] && !releasing)
            {
                activateThread(threadi);
                wait();
                waited = true;
            }
            if (releasing && waited)
            {
                break;
            }
        }

        if (lastArrival)
        {
            releasing = true;
            yieldToMaster();
        }

        if (thisThreadNum_ == 0)
        {
            threadsArrived = false;
        }
        if (thisThreadNum_ == numActiveThreads_-1)
        {
            releasing = false;
        }
    }

    // Wake up next before terminating
    activateNext();
}


// * * * * * * * * * * * * * Static functions  * * * * * * * * * * * * * * * //

void Foam::serialThreads::create(const label maxNumThreads)
{
    activeThread_ = -1;
    numThreads_ = maxNumThreads;
    numActiveThreads_ = maxNumThreads;
    thisThreadNum_ = -1;
    paused_ = 0;
    quit_ = false;

    // Create threads
    threads_ = std::vector<std::thread>(numThreads_);
    threadNotify_ = std::vector<std::condition_variable>(numThreads_);
    for (label i = 0; i < numThreads_; i++)
    {
        {
            std::lock_guard<std::mutex> modifyLock(activeThreadMutex_);
            activeThread_ = i;
        }
        threads_[i] = std::thread(threadEntryPoint, i);

        // Spawn the threads sequentially - we (the master) wait here until
        // the thread just created is ready, and only then continue to create
        // the next one (seems necessary on mingw with win32 threads)
        wait();
    }
    numActiveThreads_ = 0;
}

bool Foam::serialThreads::active()
{
    return numActiveThreads_ > 1;
}

void Foam::serialThreads::end()
{
    if (thisThreadNum_ == -1)
    {
        quit_ = true;
        numActiveThreads_ = numThreads_;

        // Wait for all threads to finish
        for (label i = 0; i < numThreads_; i++)
        {
            activateThread(i);
            threads_[i].join();
        }

        activeThread_ = -1;
        numThreads_ = 0;
        numActiveThreads_ = 0;
        thisThreadNum_ = -1;
        paused_ = 0;
    }
}

void Foam::serialThreads::yieldToNext()
{
    // Only do anything if we are running threads
    if (thisThreadNum_ != -1 && numActiveThreads_ > 1 && !paused_)
    {
        activateNext();
        wait();
    }
}

void Foam::serialThreads::yieldToMaster()
{
    if (thisThreadNum_ != -1)
    {
        activateMaster();
        wait();
    }
}

void Foam::serialThreads::level()
{
    static List<bool> threadsArrived(numThreads_, false);
    static bool releasing = false;

    // Only do anything if we are running threads
    if (numActiveThreads_ > 1 && !paused_)
    {
        threadsArrived[thisThreadNum_] = true;
        bool lastArrival = true;
        for (label threadi = 0; threadi < numActiveThreads_; threadi++)
        {
            lastArrival &= threadsArrived[threadi];
            bool waited = false;
            while (!threadsArrived[threadi] && !releasing)
            {
                activateThread(threadi);
                wait();
                waited = true;
            }
            if (releasing && waited)
            {
                break;
            }
        }

        if (lastArrival)
        {
            releasing = true;
            // Release them in order
            activateThread(0);
            wait();
        }

        if (thisThreadNum_ == 0)
        {
            threadsArrived = false;
        }
        if (thisThreadNum_ == numActiveThreads_-1)
        {
            releasing = false;
        }
    }
}


Foam::label Foam::serialThreads::thisThreadNum()
{
    if (thisThreadNum_ < 0)
    {
        FatalErrorInFunction
            << "Thread number requested outside threaded region"
            << nl << endl << exit(FatalError);
    }
    return thisThreadNum_;
}


Foam::label Foam::serialThreads::numThreads()
{
    if (thisThreadNum_ < 0)
    {
        FatalErrorInFunction
            << "Number of threads requested outside threaded region"
            << nl << endl << exit(FatalError);
    }
    return numActiveThreads_;
}


void Foam::serialThreads::pauseSwitching()
{
    paused_++;
}


void Foam::serialThreads::resumeSwitching()
{
    paused_--;
    if (paused_ < 0)
    {
        FatalErrorInFunction
            << "Switching cannot be resumed more times than paused."
            << exit(FatalError);
    }
}


// ************************************************************************* //
