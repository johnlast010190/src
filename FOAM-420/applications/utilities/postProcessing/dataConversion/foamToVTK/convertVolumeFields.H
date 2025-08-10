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
    Code chunk for converting volume and dimensioned fields
    included by foamToVTK.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{
    using reportFields = foamToVtkReportFields;

    const hashedWordList volFieldTypesWithoutInternal
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName,
    };
    const hashedWordList volFieldTypesWithInternal
    {
        volScalarField::Internal::typeName,
        volVectorField::Internal::typeName,
        volSphericalTensorField::Internal::typeName,
        volSymmTensorField::Internal::typeName,
        volTensorField::Internal::typeName
    };

    HashTable<wordHashSet> objectClasses = objects.classes();

    label nVolFields = 0;
    label nDimFields = 0;
    if (doInternal || doBoundary)
    {
        for (auto fieldType : volFieldTypesWithoutInternal)
        {
            wordHashSet& fieldNames = objectClasses(fieldType);
            label sz = fieldNames.size();
            nVolFields += sz;
        }
        for (auto fieldType : volFieldTypesWithInternal)
        {
            wordHashSet& fieldNames = objectClasses(fieldType);
            label sz = fieldNames.size();
            nVolFields += sz;
            nDimFields += sz;
        }
    }
    const hashedWordList pointFieldTypes
    {
        "pointScalarField",
        "pointVectorField",
        "pointSphericalTensorField",
        "pointSymmTensorField",
        "pointTensorField"
    };

    const label nPointFields = 0;
    if (!noPointValues)
    {
        for (auto fieldType : pointFieldTypes)
        {
            wordHashSet& fieldNames = objectClasses(fieldType);
            label sz = fieldNames.size();
            nVolFields += sz;
            nDimFields += sz;
        }
    }

    reportFields::volume(Info, objects);
    reportFields::internal(Info, objects);


    // Setup for the vtm writer.
    // For legacy format, the information added is simply ignored.

    fileName vtmOutputBase
    (
        outputDir/regionPrefix/vtkName + timeDesc
    );

    // Combined internal + boundary in a vtm file
    vtk::vtmWriter vtmWriter;

    // Collect individual boundaries into a vtm file
    vtk::vtmWriter vtmBoundaries;

    // Setup the internal writer
    autoPtr<vtk::internalWriter> internalWriter;

    // Interpolator for volume and dimensioned fields
    autoPtr<volPointInterpolation> pInterp;

    if (doInternal)
    {
        if (!noPointValues)
        {
            pInterp.reset(new volPointInterpolation(mesh));
        }

        if (vtuMeshCells.empty())
        {
            // Use the appropriate mesh (baseMesh or subMesh)
            vtuMeshCells.reset(meshProxy.mesh());
        }

        internalWriter = autoPtr<vtk::internalWriter>::New
        (
            meshProxy.mesh(),
            vtuMeshCells,
            writeOpts,
            // The output base name for internal
            (
                writeOpts.legacy()
              ? vtmOutputBase
              : (vtmOutputBase / "internal")
            ),
            Pstream::parRun()
        );

        // No sub-block for internal
        vtmWriter.append_vtu
        (
            "internal",
            vtmOutputBase.name()/"internal"
        );

        Info<< "    Internal  : "
            << internalWriter->output() << nl;

        internalWriter->writeTimeValue(mesh.time().value());
        internalWriter->writeGeometry();
    }


    // Setup the patch writers

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    PtrList<vtk::patchWriter> patchWriters;
    PtrList<PrimitivePatchInterpolation<primitivePatch>> patchInterps;

    labelList patchIds;
    if (doBoundary)
    {
        patchIds = getSelectedPatches(patches, includePatches, excludePatches);
    }

    if (oneBoundary && patchIds.size())
    {
        auto writer = autoPtr<vtk::patchWriter>::New
        (
            meshProxy.mesh(),
            patchIds,
            writeOpts,
            nearCellValue,
            // Output one patch: "boundary"
            (
                writeOpts.legacy()
              ?
                (
                    outputDir/regionPrefix/"boundary"
                  / (meshProxy.useSubMesh() ? meshProxy.name() : "boundary")
                  + timeDesc
                )
              : (vtmOutputBase / "boundary")
            ),
            Pstream::parRun()
        );

        // No sub-block for one-patch
        vtmWriter.append_vtp
        (
            "boundary",
            vtmOutputBase.name()/"boundary"
        );

        Info<< "    Boundaries: "
            << writer->output() << nl;

        writer->writeTimeValue(timeValue);
        writer->writeGeometry();

        // Transfer writer to list for later use
        patchWriters.resize(1);
        patchWriters.set(0, writer);

        // Avoid patchInterpolation for each sub-patch
        patchInterps.resize(1); // == nullptr
    }
    else if (patchIds.size())
    {
        patchWriters.resize(patchIds.size());
        if (!noPointValues)
        {
            patchInterps.resize(patchIds.size());
        }

        label nPatchWriters = 0;
        label nPatchInterps = 0;

        for (const label patchId : patchIds)
        {
            const directPolyPatch& pp = static_cast<const directPolyPatch&>(patches[patchId]);

            auto writer = autoPtr<vtk::patchWriter>::New
            (
                meshProxy.mesh(),
                labelList(one(), pp.index()),
                writeOpts,
                nearCellValue,
                // Output patch: "boundary"/name
                (
                    writeOpts.legacy()
                  ?
                    (
                        outputDir/regionPrefix/pp.name()
                      / (meshProxy.useSubMesh() ? meshProxy.name() : pp.name())
                      + timeDesc
                    )
                  : (vtmOutputBase / "boundary" / pp.name())
                ),
                Pstream::parRun()
            );

            if (!nPatchWriters)
            {
                vtmWriter.beginBlock("boundary");
                vtmBoundaries.beginBlock("boundary");
            }

            vtmWriter.append_vtp
            (
                pp.name(),
                vtmOutputBase.name()/"boundary"/pp.name()
            );

            vtmBoundaries.append_vtp
            (
                pp.name(),
                "boundary"/pp.name()
            );

            Info<< "    Boundary  : "
                << writer->output() << nl;

            writer->writeTimeValue(timeValue);
            writer->writeGeometry();

            // Transfer writer to list for later use
            patchWriters.set(nPatchWriters++, writer);

            if (patchInterps.size())
            {
                patchInterps.set
                (
                    nPatchInterps++,
                    new PrimitivePatchInterpolation<primitivePatch>(pp)
                );
            }
        }

        if (nPatchWriters)
        {
            vtmWriter.endBlock("boundary");
        }

        patchWriters.resize(nPatchWriters);
        patchInterps.resize(nPatchInterps);
    }


    // CellData
    {
        // Begin CellData
        if (internalWriter.valid())
        {
            // Optionally with cellID and procID fields
            internalWriter->beginCellData
            (
                (withMeshIds ? 1 + (internalWriter->parallel() ? 1 : 0) : 0)
              + nVolFields + nDimFields
            );

            if (withMeshIds)
            {
                internalWriter->writeCellIDs();
                internalWriter->writeProcIDs(); // parallel only
            }
        }

        if (nVolFields)
        {
            for (vtk::patchWriter& writer : patchWriters)
            {
                // Optionally with patchID field
                writer.beginCellData
                (
                    (withMeshIds ? 1 : 0)
                  + nVolFields
                );

                if (withMeshIds)
                {
                    writer.writePatchIDs();
                }
            }
        }

        writeAllVolFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true  // syncPar
        );

        writeAllDimFields
        (
            internalWriter,
            meshProxy,
            objects,
            true  // syncPar
        );

        // End CellData is implicit
    }


    // PointData
    // - only construct pointMesh on request since it constructs
    //   edge addressing
    if (!noPointValues)
    {
        // Begin PointData
        if (internalWriter.valid())
        {
            internalWriter->beginPointData
            (
                nVolFields + nDimFields + nPointFields
              + (withPointIds ? 1 : 0)
            );

            if (withPointIds)
            {
                internalWriter->writePointIDs();
            }
        }

        forAll(patchWriters, writeri)
        {
            const label nPatchFields =
            (
                (
                    writeri < patchInterps.size() && patchInterps.set(writeri)
                  ? nVolFields
                  : 0
                )
              + nPointFields
            );

            if (nPatchFields)
            {
                patchWriters[writeri].beginPointData(nPatchFields);
            }
        }

        writeAllVolFields
        (
            internalWriter, pInterp,
            patchWriters,   patchInterps,
            meshProxy,
            objects,
            true  // syncPar
        );

        writeAllDimFields
        (
            internalWriter, pInterp,
            meshProxy,
            objects,
            true  // syncPar
        );

        writeAllPointFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true  // syncPar
        );

        // End PointData is implicit
    }


    // Finish writers
    if (internalWriter.valid())
    {
        internalWriter->close();
    }

    for (vtk::patchWriter& writer : patchWriters)
    {
        writer.close();
    }

    pInterp.clear();
    patchWriters.clear();
    patchInterps.clear();


    // Collective output

    if (Pstream::master())
    {
        // Naming for vtm, file series etc.
        fileName outputName(vtmOutputBase);

        if (writeOpts.legacy())
        {
            if (doInternal)
            {
                // Add to file-series and emit as JSON

                outputName.ext(vtk::legacy::fileExtension);

                fileName seriesName(vtk::seriesWriter::base(outputName));

                vtk::seriesWriter& series = vtkSeries(seriesName);

                // First time?
                // Load from file, verify against filesystem,
                // prune time >= currentTime
                if (series.empty())
                {
                    series.load(seriesName, true, timeValue);
                }

                series.append(timeValue, timeDesc);
                series.write(seriesName);
            }
        }
        else
        {
            if (vtmWriter.size())
            {
                // Emit ".vtm"

                outputName.ext(vtmWriter.ext());

                fileName seriesName(vtk::seriesWriter::base(outputName));

                vtmWriter.setTime(timeValue);
                vtmWriter.write(outputName);

                // Add to file-series and emit as JSON

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

                // Add to multi-region vtm
                vtmMultiRegion.add(regionName, regionPrefix, vtmWriter);
            }

            if (vtmBoundaries.size())
            {
                // Emit "boundary.vtm" with collection of boundaries

                // Naming for vtm
                fileName outputName(vtmOutputBase / "boundary");
                outputName.ext(vtmBoundaries.ext());

                vtmBoundaries.setTime(timeValue);
                vtmBoundaries.write(outputName);
            }
        }
    }
}


// ************************************************************************* //
