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
    (c) 2010-2011, Esi Ltd

Application
    foamCheckMesh

Description
    Quality checks for foamHexMesh

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/Time/Time.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"
#include "regionProperties/regionProperties.H"
#include "motionSmoother/polyMeshGeometry/meshStatistics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );

    argList::validOptions.insert
    (
        "writeMetrics",
        "List of metrics: nonOrthogonality pyramids skewness weights volumeRatio determinant"
    );
    argList::validOptions.insert("writeAllMetrics", "");

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool allRegions = args.optionFound("allRegions");

    wordList regionNames;
    wordList regionDirs;
    if (allRegions)
    {
        Info<< "Reconstructing for all regions in regionProperties" << nl
            << endl;
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
        regionDirs = regionNames;
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            regionDirs = regionNames;
        }
        else
        {
            regionNames = wordList(1, fvMesh::defaultRegion);
            regionDirs = wordList(1, word::null);
        }
   }

   forAll(regionNames, regionI)
   {
        const word& regionName = regionNames[regionI];

        Info<< "\n\nCalculating quality for mesh " << regionName << nl
            << endl;

        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );

        word dictName("foamHexMeshDict");

        IOobject meshHeader
        (
            dictName,
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        );

        if (!meshHeader.typeHeaderOk<IOdictionary>(true))
        {
            dictName = "snappyHexMeshDict";
        }

        IOdictionary foamDict
        (
            IOobject
            (
                dictName,
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        const dictionary& qualityDict =
            foamDict.subDict("meshQualityControls");

        Info<< "Using mesh quality settings from meshQualityControls:" << nl
            << qualityDict << nl << endl;

        wordList metrics;

        if (args.optionFound("writeAllMetrics"))
        {
            metrics.setSize(6);
            metrics[0] = "nonOrthogonality";
            metrics[1] = "pyramids";
            metrics[2] = "skewness" ;
            metrics[3] = "weights";
            metrics[4] = "volumeRatio";
            metrics[5] = "determinant";
            Info<<"Outputting All metric fields: "<< metrics <<endl;
        }
        else if (args.optionFound("writeMetrics"))
        {
            args.optionLookup("writeMetrics")() >> metrics;
            Info<<"Outputting selected metric fields: "<< metrics <<endl;
        }

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);
            Info<< "Time = " << runTime.timeName() << endl;

            if (!timeI || mesh.readUpdate() != polyMesh::UNCHANGED)
            {
                PtrList<meshStatistics> metricsFields(0);
                if (metrics.size())
                {
                    metricsFields.setSize(metrics.size());
                    forAll(metrics, j )
                    {
                        word metricName = metrics[j];
                        dictionary histSettings;
                        if (qualityDict.found("histograms"))
                        {
                            dictionary allHist =
                                qualityDict.subDict("histograms");
                            if (allHist.found(metricName))
                            {
                                histSettings = allHist.subDict(metricName);
                            }
                        }

                        metricsFields.set
                        (
                            j,
                            new meshStatistics
                            (
                                mesh,
                                metricName,
                                histSettings
                             )
                         );
                    }
                }

                faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/1000);
                // Check with report=true
                label nErrorFaces = motionSmoother::checkMesh
                (
                    false,
                    mesh,
                    qualityDict,
                    errorFaces,
                    true,
                    &metricsFields
                );


                forAll(metrics, j)
                {
                    Info<<"Writing metric quality field: "<<metrics[j]<<endl;
                    metricsFields[j].field().write();
                    metricsFields[j].writeStats();
                }

                if (nErrorFaces > 0)
                {
                    Info<< "Writing "<< nErrorFaces
                        << " error faces to faceSet "
                        << errorFaces.objectPath() << nl << endl;
                    errorFaces.write();
                }
            }
        }
   }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
