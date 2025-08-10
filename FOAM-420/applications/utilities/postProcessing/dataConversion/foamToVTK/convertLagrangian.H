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
    (c) 2018 OpenCFD Ltd.

Description
    Code chunk for post-processing conversion of cloud(s) to VTK PolyData
    (.vtp extension).

\*---------------------------------------------------------------------------*/

if (doLagrangian)
{
    WarningInFunction
        << "Lagrangian conversion is not supported in this version of foamToVTK"
        << endl;
    /*
    const fileName cloudPrefix = (regionPrefix/cloud::prefix);

    fileNameList temporary_file_name_list = readDir
    (
        mesh.time().path()/mesh.time().timeName()/cloudPrefix,
        fileName::DIRECTORY
    );

    List<word> cloudNames(temporary_file_name_list.size());

    for (label i = 0; i < temporary_file_name_list.size(); ++i)
    {
        cloudNames[i] = temporary_file_name_list[i].name();
    }

    // Synchronise cloud names
    Pstream::combineGather(cloudNames, ListUniqueEqOp<word>());
    Pstream::combineScatter(cloudNames);

    // Consistent order
    Foam::sort(cloudNames);


    for (const word& cloudName : cloudNames)
    {
        IOobjectList cloudObjs(mesh, runTime.timeName(), cloudPrefix/cloudName);

        bool isCloud = false;
        if (cloudObjs.erase("coordinates"))
        {
            isCloud = true;
        }
        if (cloudObjs.erase("positions"))
        {
            isCloud = true;
        }

        if (!returnReduce(isCloud, orOp<bool>()))
        {
            continue;
        }

        // Limited to basic IOField types
        const hashedWordList basicFields
        {
            "labelField",               //< labelIOField
            "scalarField",              //< scalarIOField
            "vectorField",              //< vectorOField
            "sphericalTensorField",     //< sphericalTensorIOField
            "symmTensorField",          //< symmTensorIOField
            "tensorField"               //< tensorIOField
        };
        for (auto fieldType : basicFields)
        {
            cloudObjs.filterKeys
            (
                fieldType,
                true  // prune
            );
        }

        // Are there cloud fields (globally)?
        if (returnReduce(cloudObjs.empty(), andOp<bool>()))
        {
            continue;
        }

        vtk::lagrangianWriter writer
        (
            meshProxy.baseMesh(),
            cloudName,
            writeOpts,
            // Output name for the cloud
            (
                outputDir/regionPrefix/cloud::prefix
              / cloudName/cloudName + timeDesc
            ),
            Pstream::parRun()
        );

        Info<< "    Lagrangian: "
            << writer.output() << nl;

        writer.writeTimeValue(mesh.time().value());
        writer.writeGeometry();

        // Begin CellData/PointData
        writer.beginParcelData();

        writer.writeFields<label>(cloudObjs);
        writer.writeFields<scalar>(cloudObjs);
        writer.writeFields<vector>(cloudObjs);
        writer.writeFields<sphericalTensor>(cloudObjs);
        writer.writeFields<symmTensor>(cloudObjs);
        writer.writeFields<tensor>(cloudObjs);

        fileName outputName(writer.output());

        writer.close();

        if (Pstream::master())
        {
            // Add to file-series and emit as JSON

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }
    }
    */
}


// ************************************************************************* //
