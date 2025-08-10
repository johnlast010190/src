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
    (c) 2017 OpenFOAM Foundation
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "preciceAdapterSolverObject.H"

// OpenFOAM header files
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(preciceAdapterSolverObject, 0);
}
}

makeFvSolverOption(preciceAdapterSolverObject);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::preciceAdapterSolverObject::preciceAdapterSolverObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    adapter_(mesh_.time(), mesh_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::preciceAdapterSolverObject::~preciceAdapterSolverObject()
{
#ifdef ADAPTER_ENABLE_TIMINGS
    Info<< "-------------------- preCICE adapter timers (primary rank) --------------------------" << nl;
    Info<< "Total time in adapter + preCICE: " << timeInAll_.str() << " (format: day-hh:mm:ss.ms)" << nl;
    Info<< "  For setting up (S):            " << timeInSetup_.str() << " (read() function)" << nl;
    Info<< "  For all iterations (I):        " << timeInExecute_.str() << " (execute() and adjustTimeStep() functions)" << nl << nl;
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::preciceAdapterSolverObject::read(const dictionary& dict)
{}


bool Foam::fv::preciceAdapterSolverObject::initialise()
{
#ifdef ADAPTER_ENABLE_TIMINGS
    // Save the current wall clock time stamp to the clock
    clockValue clock;
    clock.update();
#endif
    Info<< "Calling configure " << endl;
    adapter_.configure();

#ifdef ADAPTER_ENABLE_TIMINGS
    // Accumulate the time in this section into a global timer.
    // Same in all function object methods.
    timeInAll_ += clock.elapsed();
    timeInSetup_ = clock.elapsed();
#endif

    return true;
}


void Foam::fv::preciceAdapterSolverObject::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("readPreciceData");
    solveNames.append("preciceSolve");
    correctorMembers.insert
    (
        "preciceCouplingCorrector",
        {"readPreciceData", "preciceSolve", solverObject::outerCorrectorName}
    );

    requiredDependencies.insert
    (
        "fvMesh",
        {"readPreciceData"}
    );
    optionalDependencies.insert
    (
        "preciceSolve",
        {"p"}
    );
}


bool Foam::fv::preciceAdapterSolverObject::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == "preciceCouplingCorrector")
    {
        return
            !adapter_.isCouplingOngoing()
         || (corrector != 0 && adapter_.isCouplingTimeWindowComplete());
    }
    else
    {
        return true;
    }
}

void Foam::fv::preciceAdapterSolverObject::correct
(
    const word& solveName,
    const word& regionName
)
{
    // preCICE has already moved on to the next time step by the time we know
    // it has converged. Therefore, don't call preCICE in the final iteration
    if (!finalIter_["preciceCouplingCorrector"] && adapter_.isCouplingOngoing())
    {
        if (solveName == "readPreciceData")
        {
#ifdef ADAPTER_ENABLE_TIMINGS
            clockValue clock;
            clock.update();
#endif

            adapter_.readCouplingData();

#ifdef ADAPTER_ENABLE_TIMINGS
            timeInAll_ += clock.elapsed();
#endif
        }
        else
        {
#ifdef ADAPTER_ENABLE_TIMINGS
            clockValue clock;
            clock.update();
#endif

            adapter_.execute();

#ifdef ADAPTER_ENABLE_TIMINGS
            timeInAll_ += clock.elapsed();
            timeInExecute_ += clock.elapsed();
    #endif
        }
    }
}


void Foam::fv::preciceAdapterSolverObject::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if
    (
        correctorName == solverObject::timeLoopCorrectorName
        && finalIter
    )
    {
        // Call finalize here in case the simulation completed before the
        // coupling completed
        adapter_.finalize();
    }
}

// ************************************************************************* //
