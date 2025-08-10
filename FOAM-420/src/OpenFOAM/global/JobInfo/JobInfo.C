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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "global/JobInfo/JobInfo.H"
#include "include/OSspecific.H"
#include "global/clock/clock.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "global/foamVersion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::JobInfo::writeJobInfo(Foam::debug::infoSwitch("writeJobInfo", 0));
Foam::JobInfo Foam::jobInfo;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::JobInfo::write(Ostream& os) const
{
    if (writeJobInfo && Pstream::master())
    {
        if (os.good())
        {
            dictionary::write(os, false);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


void Foam::JobInfo::end(const word& terminationType)
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        add("cpuTime", cpuTime_.elapsedCpuTime());
        add("endDate", clock::date());
        add("endTime", clock::clockTime());

        if (!found("termination"))
        {
            add("termination", terminationType);
        }

        Foam::rm(runningDir_/jobFileName_);
        write(OFstream(finishedDir_/jobFileName_)());
    }

    constructed = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JobInfo::JobInfo()
:
    jobFileName_(),
    runningDir_(),
    finishedDir_(),
    cpuTime_()
{
    name() = "JobInfo";

    if (writeJobInfo && Pstream::master())
    {
        string jobDir = getEnv("FOAM_JOB_DIR");
        if (jobDir.empty())
        {
            // Fallback: ~/.FOAMcore/jobControl
            jobDir = home()/FOAM_USER_RESOURCE_DIRNAME/"jobControl";
        }

        jobFileName_ = hostName() + '.' + Foam::name(pid());
        runningDir_  = jobDir/"runningJobs";
        finishedDir_ = jobDir/"finishedJobs";

        if (!isDir(jobDir) && !mkDir(jobDir))
        {
            FatalErrorInFunction
                << "No JobInfo directory:  FOAM_JOB_DIR=" << jobDir
                << Foam::exit(FatalError);
        }
        if (!isDir(runningDir_) && !mkDir(runningDir_))
        {
            FatalErrorInFunction
                << "No JobInfo directory:  " << runningDir_
                << Foam::exit(FatalError);
        }
        if (!isDir(finishedDir_) && !mkDir(finishedDir_))
        {
            FatalErrorInFunction
                << "No JobInfo directory:  " << finishedDir_
                << Foam::exit(FatalError);
        }
    }

    constructed = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::JobInfo::~JobInfo()
{
    signalEnd();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JobInfo::write() const
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        const fileName output = runningDir_/jobFileName_;
        if (!write(OFstream(output)()))
        {
            FatalErrorInFunction
                << "Failed to write to JobInfo file " << output
                << Foam::exit(FatalError);
        }
    }
}


void Foam::JobInfo::end()
{
    end("normal");
}


void Foam::JobInfo::exit()
{
    end("exit");
}


void Foam::JobInfo::abort()
{
    end("abort");
}


void Foam::JobInfo::signalEnd() const
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        Foam::mv(runningDir_/jobFileName_, finishedDir_/jobFileName_);
    }
    constructed = false;
}


// ************************************************************************* //
