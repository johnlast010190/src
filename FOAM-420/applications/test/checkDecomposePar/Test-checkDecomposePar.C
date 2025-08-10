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
    (c) 2017-2021 OpenCFD Ltd.

Description
    Check decomposition from kaffpa (KaHIP) output.
    foamToMetisGraph was likely used for producing the kaffpa input.

\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cpuTime/cpuTime.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "regionProperties/regionProperties.H"
#include "decompositionInformation.H"
#include "decompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check decomposition from kaffpa (KaHIP) output"
    );

    argList::noParallel();
    argList::noBanner();

    argList::addOption
    (
        "decomposeParDict",
        "file",
        "read decomposePar dictionary from specified location"
    );

    Foam::argList::addBoolOption
    (
        "allRegions",
        "Use all regions in regionProperties"
    );

    Foam::argList::addOption
    (
        "region",
        "name",
        "Use specified mesh region. Eg, -region gas"
    );

    argList::addBoolOption
    (
        "verbose",
        "more information about decomposition"
    );

    argList::validArgs.append("kaffpa-output-file");

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    const fileName decompFile(args[1]);
    const bool verbose = args.optionFound("verbose");

    // Set time from database
    #include "createTime.H"

    // Allow override of time
    instantList times = timeSelector::selectIfPresent(runTime, args);

    // Allow override of decomposeParDict location
    fileName decompDictFile;
    args.optionReadIfPresent("decomposeParDict", decompDictFile);

    // Get region names
    //#include "getAllRegionOptions.H"
    wordList regionNames;
    {
        wordRes selectByName;

        if (args.optionFound("allRegions"))
        {
            regionProperties rp(runTime);
            forAllConstIter(HashTable<wordList>, rp, iter)
            {
                const wordList& regions = iter();
                forAll(regions, i)
                {
                    if (findIndex(regionNames, regions[i]) == -1)
                    {
                        regionNames.append(regions[i]);
                    }
                }
            }

            if (regionNames.empty())
            {
                Info<< "Warning: No regionProperties, assume default region"
                    << nl << endl;
            }
            else
            {
                Info<< "Using all regions: " << flatOutput(regionNames) << nl;
            }
        }
        else
        {
            word regionName;
            if (args.optionReadIfPresent("region", regionName))
            {
                regionNames = wordList(1, regionName);
            }
            else
            {
                regionNames = wordList(1, polyMesh::defaultRegion);
            }
        }
    }
    labelList cellToProc;

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        // const word& regionDir =
        // (
        //     regionName != polyMesh::defaultRegion
        //   ? regionName
        //   : word::null
        // );

        Info<< "\n\nDecomposing mesh " << regionName << nl << endl;
        Info<< "Create mesh..." << flush;

        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Info<< " nCells = " << mesh.nCells() << endl;

        // Expected format is a simple ASCII list
        cellToProc.setSize(mesh.nCells());
        {
            IFstream is(decompFile);

            forAll(cellToProc, celli)
            {
                cellToProc[celli] = readLabel(is);
            }
        }

        const label nDomains = max(cellToProc) + 1;

        CompactListList<label> cellCells;
        decompositionMethod::calcCellCells
        (
            mesh,
            identity(mesh.nCells()),
            mesh.nCells(),
            false,
            cellCells
        );

        decompositionInformation info
        (
            cellCells,
            cellToProc,
            nDomains
        );

        if (verbose)
        {
            info.printDetails(Info);
            Info<< nl;
        }
        info.printSummary(Info);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
