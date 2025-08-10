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
    (c) 2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "processorRunTimes.H"
#include "decompositionModel.H"
#include "db/Time/timeSelector.H"
#include "global/argList/argList.H"
#include "regionProperties/regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorRunTimes::processorRunTimes
(
    const word& name,
    const argList& args
)
:
    completeRunTime_(name, args),
    decomposeDictFile_(),
    nProcs_
    (
        fileHandler().nProcs
        (
            args.path(),
            regionDir
            (
                selectRegionNames
                (
                    args,
                    Time
                    (
                        Time::controlDictName,
                        completeRunTime_.rootPath(),
                        completeRunTime_.caseName()
                       /fileName(word("processor0"))
                    )
                )[0]
            )
        )
    ),
    procRunTimes_(nProcs_)
{
    forAll(procRunTimes_, proci)
    {
        procRunTimes_.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                completeRunTime_.rootPath(),
                completeRunTime_.caseName()
               /fileName(word("processor") + Foam::name(proci))
            )
        );

        procRunTimes_[proci].setTime(completeRunTime_);
    }
}


Foam::processorRunTimes::processorRunTimes
(
    const word& name,
    const argList& args,
    const fileName& decomposeDictFile
)
:
    completeRunTime_(name, args),
    decomposeDictFile_(decomposeDictFile),
    nProcs_
    (
        readLabel
        (
            IOdictionary
            (
                decompositionModel::selectIO
                (
                    IOobject
                    (
                        "decomposeParDict",
                        completeTime().time().system(),
                        completeTime(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    decomposeDictFile
                )
            ).lookup("numberOfSubdomains")
        )
    ),
    procRunTimes_(nProcs_)
{
    forAll(procRunTimes_, proci)
    {
        procRunTimes_.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                completeRunTime_.rootPath(),
                completeRunTime_.caseName()
               /fileName(word("processor") + Foam::name(proci))
            )
        );

        procRunTimes_[proci].setTime(completeRunTime_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorRunTimes::~processorRunTimes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::processorRunTimes::setTime
(
    const instant& inst,
    const label newIndex
)
{
    completeRunTime_.setTime(inst, newIndex);

    forAll(procRunTimes_, proci)
    {
        procRunTimes_[proci].setTime(inst, newIndex);
    }
}


Foam::instantList Foam::processorRunTimes::selectComplete(const argList& args)
{
    instantList timeDirs =
        timeSelector::selectIfPresent(completeRunTime_, args);

    forAll(procRunTimes_, proci)
    {
        procRunTimes_[proci].setTime(completeRunTime_);
    }

    return timeDirs;
}


Foam::instantList Foam::processorRunTimes::selectProc(const argList& args)
{
    instantList timeDirs =
        timeSelector::select0(procRunTimes_[0], args);

    completeRunTime_.setTime(procRunTimes_[0]);

    for (label proci = 1; proci < nProcs(); proci ++)
    {
        procRunTimes_[proci].setTime(procRunTimes_[0]);
    }

    return timeDirs;
}


// ************************************************************************* //
