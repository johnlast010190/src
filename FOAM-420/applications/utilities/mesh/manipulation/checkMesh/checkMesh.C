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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.

Application
    checkMesh

Group
    grpMeshManipulationUtilities

Description
    Checks validity of a mesh.

Usage
    \b checkMesh [OPTION]

    Options:
      - \par -allGeometry
        Checks all (including non finite-volume specific) geometry

      - \par -allTopology
        Checks all (including non finite-volume specific) addressing

      - \par -meshQuality
        Checks against user defined (in \a system/meshQualityDict) quality
        settings

      - \par -region \<name\>
        Specify an alternative mesh region.

    \param -writeSets \<surfaceFormat\> \n
    Reconstruct all cellSets and faceSets geometry and write to postProcessing
    directory according to surfaceFormat (e.g. vtk or ensight). Additionally
    reconstructs all pointSets and writes as vtk format.

    \param -writeAllFields \n
    Writes all additional mesh fields.

    \param -writeFields '(\<fieldName\>)' \n
    Writes selected mesh fields.

    \param -writeAllMetrics \n
    Writes all mesh quality fields.

    \param -writeMetrics '(\<fieldName\>)' \n
    Writes all mesh quality fields.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "sampledSurface/writers/surfaceWriter.H"
#include "sampledSetWriters/vtk/vtkSetWriter.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "regionProperties/regionProperties.H"
#include "meshRefinement/meshRefinement.H"

#include "checkTools.H"
#include "checkTopology.H"
#include "checkGeometry.H"
#include "checkMeshQuality.H"
#include "writeFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "include/addRegionOption.H"
    argList::addOption
    (
        "regions",
        "(name1 .. nameN)",
        "specify list of mesh regions"
    );
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "noTopology",
        "skip checking the mesh topology"
    );
    argList::addBoolOption
    (
        "noGeometry",
        "skip checking the mesh geometry"
    );
    argList::addBoolOption
    (
        "allGeometry",
        "include bounding box checks"
    );
    argList::addBoolOption
    (
        "allTopology",
        "include extra topology checks"
    );
    argList::addBoolOption
    (
        "writeAllFields",
        "write all additional mesh volFields : cellShapes, cellVolume, cellCentres"
    );
    argList::addOption
    (
        "writeFields",
        "wordList",
        "write list additional mesh fields : cellShapes, cellVolume, cellCentres"
    );
    argList::addBoolOption
    (
        "meshQuality",
        "read user-defined mesh quality criterions from system/meshQualityDict"
    );
    argList::addOption
    (
        "writeSets",
        "surfaceFormat",
        "reconstruct and write all faceSets and cellSets in selected format"
    );
    argList::addOption
    (
        "writeMetrics",
        "List of metrics: nonOrthogonality pyramids skewness weights volumeRatio determinant warpage aspectRatio"
    );
    argList::addBoolOption
    (
        "writeAllMetrics",
        "Write fields : nonOrthogonality pyramids skewness weights volumeRatio determinant warpage aspectRatio"
    );
    argList::addBoolOption
    (
        "writeCombinedMetric",
        "Write combined metric field of all fields specified by writeMetrics or writeAllMetrics"
    );
    argList::addOption
    (
        "sourceTargetBoundaryMatching",
        "List<Tuple2<wordReList, wordReList>>",
        "Calculate boundary matching error - eg '( ((amiIn) (amiOut)) )' "
    );
    argList::addBoolOption
    (
        "writeBoundaryError",
        "Write source/target boundary error vtk file"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool noTopology  = args.optionFound("noTopology");
    const bool noGeometry  = args.optionFound("noGeometry");
    const bool allGeometry = args.optionFound("allGeometry");
    const bool allTopology = args.optionFound("allTopology");
    const bool meshQuality = args.optionFound("meshQuality");
    const bool writeAllMetrics = args.optionFound("writeAllMetrics");
    const bool writeMetrics = args.optionFound("writeMetrics");
    const bool writeCombinedMetric = args.optionFound("writeCombinedMetric");
    const bool writeErrorVTK  = args.optionFound("writeBoundaryError");
    List<Tuple2<List<wordRe>, List<wordRe>>> sourceTargetMatching;
    args.optionReadIfPresent
    (
        "sourceTargetBoundaryMatching",
        sourceTargetMatching
    );

    wordList metrics;

    if
    (
        writeAllMetrics || (writeCombinedMetric && !writeMetrics)
    )
    {
        metrics.setSize(8);
        metrics[0] = "nonOrthogonality";
        metrics[1] = "pyramids";
        metrics[2] = "skewness" ;
        metrics[3] = "weights";
        metrics[4] = "volumeRatio";
        metrics[5] = "determinant";
        metrics[6] = "warpage";
        metrics[7] = "aspectRatio";
        Info<<"Outputting All metric fields: "<< metrics <<endl;
    }
    else if (writeMetrics)
    {
        args.optionLookup("writeMetrics")() >> metrics;
        Info<<"Outputting selected metric fields: "<< metrics <<endl;
    }

    word surfaceFormat;
    const bool writeSets = args.optionReadIfPresent("writeSets", surfaceFormat);
    HashSet<word> selectedFields;
    bool writeFields = args.optionReadIfPresent
    (
        "writeFields",
        selectedFields
    );
    if (!writeFields && args.optionFound("writeAllFields"))
    {
        selectedFields.insert("cellShapes");
        selectedFields.insert("cellVolume");
        selectedFields.insert("cellCentres");
    }

    if (noTopology)
    {
        Info<< "Disabling all topology checks." << nl << endl;
    }
    if (noGeometry)
    {
        Info<< "Disabling all geometry checks." << nl << endl;
    }
    if (allTopology)
    {
        Info<< "Enabling all (cell, face, edge, point) topology checks."
            << nl << endl;
    }
    if (allGeometry)
    {
        Info<< "Enabling all geometry checks." << nl << endl;
    }
    if (meshQuality)
    {
        Info<< "Enabling user-defined geometry checks." << nl << endl;
    }
    if (writeSets)
    {
        Info<< "Reconstructing and writing " << surfaceFormat
            << " representation"
            << " of all faceSets and cellSets." << nl << endl;
    }

    wordList regionNames;
    if (args.optionFound("allRegions"))
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
    }
    else if (args.optionFound("regions"))
    {
        regionNames = args.optionRead<wordList>("regions");
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

    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];

        Info<< "Checking mesh for region : " << regionName << nl <<endl;

        Foam::fvMesh mesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );

        autoPtr<IOdictionary> qualDict;

        IOobject meshHeader
        (
            "meshQualityDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        );

        if (meshHeader.typeHeaderOk<IOdictionary>(true))
        {
            qualDict.reset
            (
                new IOdictionary
                (
                    IOobject
                    (
                        "meshQualityDict",
                        mesh.time().system(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                 )
            );
        }
        else
        {
            IOobject foamHeader
            (
                "foamHexMeshDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );

            if (foamHeader.typeHeaderOk<IOdictionary>(true))
            {
                IOdictionary foamDict
                (
                    IOobject
                    (
                        "foamHexMeshDict",
                        mesh.time().system(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                const dictionary& qualityDict =
                    foamDict.subDict("meshQualityControls");
                qualDict.reset
                (
                    new IOdictionary
                    (
                        IOobject
                        (
                            "meshQualityDict",
                            mesh.time().system(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        qualityDict
                     )
                );
            }
            else
            {
                dictionary metricDict("meshQualityControls");
                forAll(metrics, metricI)
                {
                    word metricName = metrics[metricI];
                    if (metricName == "nonOrthogonality")
                    {
                        metricDict.add("maxNonOrtho", scalar(70));
                    }
                    else if (metricName == "pyramids")
                    {
                        metricDict.add("minVol", scalar(1e-16));
                    }
                    else if (metricName == "skewness")
                    {
                        metricDict.add("maxBoundarySkewness", scalar(20));
                        metricDict.add("maxInternalSkewness", scalar(6));
                    }
                    else if (metricName == "weights")
                    {
                        metricDict.add("minFaceWeight", scalar(0.08));
                    }
                    else if (metricName == "volumeRatio")
                    {
                        metricDict.add("minVolRatio", scalar(0.02));
                    }
                    else if (metricName == "determinant")
                    {
                        metricDict.add("minDeterminant", scalar(0.0001));
                    }
                    else if (metricName == "aspectRatio")
                    {
                        metricDict.add("maxCellAspectRatio", scalar(1000));
                    }
                    else if (metricName == "warpage")
                    {
                        metricDict.add("maxInternalWarpage", scalar(0.5));
                        metricDict.add("maxBoundaryWarpage", scalar(0.5));
                    }
                }

                qualDict.reset
                (
                    new IOdictionary
                    (
                        IOobject
                        (
                            "meshQualityDict",
                            mesh.time().system(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                         ),
                        metricDict
                     )
                );
            }
        }

        //Check all output metrics are activated
        forAll(metrics, metricI)
        {
            word metricName = metrics[metricI];
            if (metricName == "nonOrthogonality")
            {
                if (!qualDict().found("maxNonOrtho"))
                {
                    qualDict().add("maxNonOrtho", scalar(70));
                }
                else
                {
                    scalar maxOrtho =
                        readScalar(qualDict().lookup("maxNonOrtho"));
                    if (maxOrtho > 180.0-SMALL)
                    {
                        qualDict().add("maxNonOrtho", scalar(70),true);
                    }
                }
            }
            else if (metricName == "pyramids")
            {
                if (!qualDict().found("minVol"))
                {
                    qualDict().add("minVol", scalar(1e-16));
                }
                else
                {
                    scalar minVol =
                        readScalar(qualDict().lookup("minVol"));
                    if (minVol == -GREAT)
                    {
                        qualDict().add("minVol",scalar(1e-16),true);
                    }
                }
            }
            else if (metricName == "skewness")
            {
                if (!qualDict().found("maxBoundarySkewness"))
                {
                    qualDict().add("maxBoundarySkewness", scalar(20));
                }
                else
                {
                    scalar maxBoundSkew =
                        readScalar(qualDict().lookup("maxBoundarySkewness"));
                    if (maxBoundSkew < 0)
                    {
                        qualDict().add("maxBoundarySkewness",scalar(20),true);
                    }
                }

                if (!qualDict().found("maxInternalSkewness"))
                {
                    qualDict().add("maxInternalSkewness", scalar(6));
                }
                else
                {
                    scalar maxIntSkew =
                        readScalar(qualDict().lookup("maxInternalSkewness"));
                    if (maxIntSkew < 0)
                    {
                        qualDict().add("maxInternalSkewness",scalar(6),true);
                    }
                }
            }
            else if (metricName == "weights")
            {
                if (!qualDict().found("minFaceWeight"))
                {
                    qualDict().add("minFaceWeight", scalar(0.08));
                }
                else
                {
                    scalar minFaceWeight =
                        readScalar(qualDict().lookup("minFaceWeight"));
                    if (minFaceWeight < 0)
                    {
                        qualDict().add("minFaceWeight",scalar(0.08),true);
                    }
                }
            }
            else if (metricName == "volumeRatio")
            {
                if (!qualDict().found("minVolRatio"))
                {
                    qualDict().add("minVolRatio", scalar(0.02));
                }
                else
                {
                    scalar minVolRatio =
                        readScalar(qualDict().lookup("minVolRatio"));
                    if (minVolRatio < 0)
                    {
                        qualDict().add("minVolRatio",scalar(0.02),true);
                    }
                }
            }
            else if (metricName == "determinant")
            {
                if (!qualDict().found("minDeterminant"))
                {
                    qualDict().add("minDeterminant", scalar(0.0001));
                }
                else
                {
                    scalar minDet =
                        readScalar(qualDict().lookup("minDeterminant"));
                    if (minDet < 0)
                    {
                        qualDict().add("minDeterminant",scalar(0.0001),true);
                    }
                }
            }
            else if (metricName == "warpage")
            {
                if (!qualDict().found("maxInternalWarpage"))
                {
                    qualDict().add("maxInternalWarpage", scalar(0.5));
                }
                else
                {
                    scalar intWarp =
                        readScalar(qualDict().lookup("maxInternalWarpage"));
                    if (intWarp <= 0)
                    {
                        qualDict().add("maxInternalWarpage",scalar(0.5),true);
                    }
                }

                if (!qualDict().found("maxBoundaryWarpage"))
                {
                    qualDict().add("maxBoundaryWarpage", scalar(0.5));
                }
                else
                {
                    scalar boundWarp =
                        readScalar(qualDict().lookup("maxBoundaryWarpage"));
                    if (boundWarp <= 0)
                    {
                        qualDict().add("maxBoundaryWarpage",scalar(0.5),true);
                    }
                }
            }
            else if (metricName == "aspectRatio")
            {
                if (!qualDict().found("maxCellAspectRatio"))
                {
                    qualDict().add("maxCellAspectRatio", scalar(1000));
                }
                else
                {
                    scalar maxAR =
                        readScalar(qualDict().lookup("maxCellAspectRatio"));
                    if (maxAR < 0)
                    {
                        qualDict().add("maxCellAspectRatio",scalar(1000),true);
                    }
                }
            }
        }

        autoPtr<surfaceWriter> surfWriter;
        autoPtr<writer<scalar>> setWriter;
        if (writeSets)
        {
            surfWriter = surfaceWriter::New(surfaceFormat);
            setWriter = writer<scalar>::New(vtkSetWriter<scalar>::typeName);
        }


        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            polyMesh::readUpdateState state = mesh.readUpdate();

            if
            (
                !timeI
                || state == polyMesh::TOPO_CHANGE
                || state == polyMesh::TOPO_PATCH_CHANGE
             )
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;

                // Reconstruct globalMeshData
                mesh.globalData();

                printMeshStats(mesh, allTopology);

                if (sourceTargetMatching.size())
                {
                    forAll(sourceTargetMatching, sti)
                    {
                        meshRefinement::calcSourceTargetMatch
                        (
                           mesh,
                           sourceTargetMatching[sti].first(),
                           sourceTargetMatching[sti].second(),
                           writeErrorVTK
                        );
                    }
                }

                label nFailedChecks = 0;

                if (!noTopology)
                {
                    nFailedChecks += checkTopology
                    (
                        mesh,
                        allTopology,
                        allGeometry,
                        surfWriter,
                        setWriter
                    );
                }

                if (!noGeometry)
                {
                    nFailedChecks += checkGeometry
                    (
                        mesh,
                        allGeometry,
                        surfWriter,
                        setWriter
                    );
                }

                if
                (
                    meshQuality || writeAllMetrics
                    || writeMetrics || writeCombinedMetric
                )
                {
                    label nFailedQualityChecks =
                        checkMeshQuality
                        (
                            mesh,
                            qualDict(),
                            metrics,
                            meshQuality,
                            writeCombinedMetric,
                            surfWriter
                        );

                    if (meshQuality)
                    {
                        nFailedChecks += nFailedQualityChecks;
                    }
                }


                // Note: no reduction in nFailedChecks necessary since is
                // counter of checks, not counter of failed cells,faces etc.

                if (nFailedChecks == 0)
                {
                    if (!noGeometry)
                    {
                        Info<< "\nMesh OK.\n" << endl;
                    }
                }
                else
                {
                    Info<< "\nWarning in " << nFailedChecks << " mesh checks.\n"
                        << endl;
                }

                // Write selected fields
                Foam::writeFields(mesh, selectedFields);
            }
            else if (state == polyMesh::POINTS_MOVED)
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;

                label nFailedChecks = 0;
                if (!noGeometry)
                {
                    nFailedChecks += checkGeometry
                    (
                        mesh,
                        allGeometry,
                        surfWriter,
                        setWriter
                    );
                }

                if
                (
                    meshQuality || writeAllMetrics
                    || writeMetrics || writeCombinedMetric
                )
                {
                    label nFailedQualityChecks = checkMeshQuality
                    (
                        mesh,
                        qualDict(),
                        metrics,
                        meshQuality,
                        writeCombinedMetric,
                        surfWriter
                    );

                    if (meshQuality)
                    {
                        nFailedChecks += nFailedQualityChecks;
                    }
                }

                if (nFailedChecks)
                {
                    if (!noGeometry)
                    {
                        Info<< "\nWarning in " << nFailedChecks
                            << " mesh checks.\n"
                            << endl;
                    }
                }
                else
                {
                    Info<< "\nMesh OK.\n" << endl;
                }

                // Write selected fields
                Foam::writeFields(mesh, selectedFields);
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
