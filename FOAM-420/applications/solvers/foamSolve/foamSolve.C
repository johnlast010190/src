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

Application
    foamSolve.C

Description
    General single- or multi-region runner for solver objects

SeeAlso
    Foam::foamSolve

\*---------------------------------------------------------------------------*/

#include "foamSolve.H"
#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add postProcessing option
    Foam::argList::addBoolOption
    (
        argList::postProcessOptionName,
        "Execute functionObjects only"
    );
    Foam::argList::addOption
    (
        "instances",
        "list",
        "specify a list of instances to run"
    );
    if (argList::postProcess(argc, argv))
    {
        Foam::timeSelector::addOptions();
        #include "include/addRegionOption.H"
        #include "include/addFunctionObjectOptions.H"

        // Set functionObject post-processing mode
        functionObject::postProcess = true;
    }

    // Add profiling option
    #if !defined( WIN32 ) && !defined( WIN64 )
    #include "include/addProfilingOption.H"
    #endif

    // Set root case
    #include "include/setRootCase.H"

    fileName exeName(args.executable());

    if (exeName == "foamCHT")
    {
        DeprecationWarningInFunction("foamCHT", "solver", 40100) << endl;
    }
    else if (exeName == "foamCoupledCompressible")
    {
        DeprecationWarningInFunction
        (
            "foamCoupledCompressible", "solver", 40100
        ) << endl;
    }

    return foamSolve().runSolver(args, exeName.nameLessExt());
}

// ************************************************************************* //
