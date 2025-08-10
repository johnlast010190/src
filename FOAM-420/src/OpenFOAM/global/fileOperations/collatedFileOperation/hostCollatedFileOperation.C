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
    (c) 2017-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "global/fileOperations/collatedFileOperation/hostCollatedFileOperation.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/Lists/PackedList/PackedBoolList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostCollatedFileOperation,
        word
    );

    // Register initialisation routine. Signals need for threaded mpi and
    // handles command line arguments
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        hostCollatedFileOperationInitialise,
        word,
        hostCollated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fileOperations::hostCollatedFileOperation::subRanks
(
    const label n
)
{
    DynamicList<label> subRanks(64);

    string ioRanksString(getEnv("FOAM_IORANKS"));
    if (!ioRanksString.empty())
    {
        IStringStream is(ioRanksString);
        labelList ioRanks(is);

        if (findIndex(ioRanks, 0) == -1)
        {
            FatalErrorInFunction
                << "Rank 0 (master) should be in the IO ranks. Currently "
                << ioRanks << exit(FatalError);
        }

        // The lowest numbered rank is the IO rank
        PackedBoolList isIOrank(n);
        isIOrank.set(ioRanks);

        for (label proci = Pstream::myProcNo(); proci >= 0; --proci)
        {
            if (isIOrank[proci])
            {
                // Found my master. Collect all processors with same master
                subRanks.append(proci);
                for
                (
                    label rank = proci+1;
                    rank < n && !isIOrank[rank];
                    ++rank
                )
                {
                    subRanks.append(rank);
                }
                break;
            }
        }
    }
    else
    {
        // Normal operation: one lowest rank per hostname is the writer
        const string myHostName(hostName());

        stringList hosts(Pstream::nProcs());
        hosts[Pstream::myProcNo()] = myHostName;
        Pstream::allGatherList(hosts);

        // Collect procs with same hostname
        forAll(hosts, proci)
        {
            if (hosts[proci] == myHostName)
            {
                subRanks.append(proci);
            }
        }
    }
    return labelList(subRanks, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    bool verbose
)
:
    collatedFileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            subRanks(Pstream::nProcs())
        ),
        (Pstream::parRun() ? labelList(0) : ioRanks()), // processor dirs
        typeName,
        verbose
    )
{
    if (verbose)
    {
        // Print a bit of information
        stringList ioRanks(Pstream::nProcs());
        if (Pstream::master(comm_))
        {
            ioRanks[Pstream::myProcNo()] = hostName()+"."+name(pid());
        }
        Pstream::gatherList(ioRanks);

        Info<< "         IO nodes:" << nl;
        for (const string& ranks : ioRanks)
        {
            if (!ranks.empty())
            {
                Info<< "             " << ranks << nl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostCollatedFileOperation::~hostCollatedFileOperation()
{
    if (comm_ != -1 && comm_ != UPstream::worldComm)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// ************************************************************************* //
