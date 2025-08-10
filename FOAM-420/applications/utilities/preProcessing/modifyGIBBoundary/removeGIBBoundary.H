    boolList unsnapPoints(mesh.points().size(), false);
    boolList checkZones(mesh.faceZones().size(), false);

    {
        // Old and new patches.
        DynamicList<polyPatch*> allPatches(patches.size());

        label startFaceI = mesh.nInternalFaces();

        labelHashSet includePatches(patches.size());
        includePatches = patches.patchSet
        (
            wordReList(args.optionLookup("names")())
        );


        // Copy old patches.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (isA<directPolyPatch>(pp))
            {
                const directPolyPatch& dpp = refCast<const
                    directPolyPatch>(pp);

                allPatches.append
                (
                    dpp.clone
                    (
                        patches,
                        patchI,
                        pp.size(),
                        startFaceI
                    ).ptr()
                );
                startFaceI += pp.size();
            }
            else
            {
                const indirectPolyPatch& ipp = refCast<const
                    indirectPolyPatch>(pp);
                if
                (
                    !includePatches[patchI]
                )
                {
                    allPatches.append
                    (
                        ipp.clone
                        (
                            patches
                        ).ptr()
                    );
                }
                else
                {

                    if (mesh.faceZones().size() != 0)
                    {
                        label fzId = ipp.zoneId();
                        word zoneName = mesh.faceZones()[fzId].name();
                        label zId = mesh.faceZones().findZoneID(zoneName);
                        if (zId != -1)
                        {
                            checkZones[zId] = true;
                        }
                    }
                    Info<< "Removing GIB boundary " << pp.name() << endl;

                    const labelList& ippPoints = ipp.meshPoints();
                    forAll(ippPoints, pI)
                    {
                        const label gpI = ippPoints[pI];
                        unsnapPoints[gpI] = true;
                    }
                }
            }
        }
        allPatches.shrink();
        mesh.removeFvBoundary();
        mesh.addFvPatches(allPatches);
    }
    patches.write();


    IOobject ioBasPoints
    (
        "basePoints",
        mesh.time().findInstance
            (
                mesh.meshDir(),
                "basePoints",
                IOobject::READ_IF_PRESENT
            ),
        Foam::fvMesh::meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (exists(ioBasPoints.objectPath()))
    {
        pointIOField basepoints
        (
            ioBasPoints
        );
        pointIOField cpoints
        (
            IOobject
            (
                "points",
                mesh.time().findInstance
                    (
                        mesh.meshDir(),
                        "points",
                        IOobject::READ_IF_PRESENT
                    ),
                Foam::fvMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        //- check and remove faceZones created on the fly.
        {
            DynamicList<word> removingFaceZones(mesh.faceZones().size());
            DynamicList<word> removingCellZones(mesh.faceZones().size());
            forAll(checkZones, zI)
            {
                if (checkZones[zI])
                {
                    const labelList& fz = mesh.faceZones()[zI];
                    bool samePoints = true;
                    forAll(fz, fI)
                    {
                        const label gfI = fz[fI];
                        forAll(mesh.faces()[gfI], pI)
                        {
                            label gpI = mesh.faces()[gfI][pI];
                            if (cpoints[gpI] != basepoints[gpI])
                            {
                                samePoints = false;
                            }
                        }
                    }
                    reduce(samePoints, andOp<bool>());
                    if (!samePoints)
                    {
                        Info<< "faceZone:" << tab
                             << mesh.faceZones()[zI].name()
                             << "removed from files"
                             << endl;
                        removingFaceZones.append(mesh.faceZones()[zI].name());
                    }
                    removingCellZones.append
                    (
                        "inactive_"+mesh.faceZones()[zI].name()
                    );
                }
            }
            removingFaceZones.shrink();
            removingCellZones.shrink();
            forAll(removingFaceZones, zI)
            {
                removeZone(mesh.faceZones(), removingFaceZones[zI]);
            }
            forAll(removingCellZones, zI)
            {
                removeZone(mesh.cellZones(), removingCellZones[zI]);
            }
        }

        //- unsnap points of deleted GIB
        forAll(unsnapPoints, pI)
        {
            if (unsnapPoints[pI])
            {
                cpoints[pI] = basepoints[pI];
            }
        }
        cpoints.write();

        IOobject meshPhiIO
        (
            "meshPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        );
        if (exists(meshPhiIO.objectPath()))
        {
            rm(meshPhiIO.filePath());
        }

        IOobject phiIO
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        );
        if (exists(phiIO.objectPath()))
        {
            rm(phiIO.filePath());
        }

        IOobject UfIO
        (
            "Uf",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        );
        if (exists(UfIO.objectPath()))
        {
            rm(UfIO.filePath());
        }


        IOobject ioHitIndex
        (
            "hitIndex",
            mesh.time().findInstance
                (
                    mesh.meshDir(),
                    "hitIndex",
                    IOobject::READ_IF_PRESENT
                 ),
             Foam::fvMesh::meshSubDir,
             mesh,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         );
        if (exists(ioHitIndex.objectPath()))
        {
            rm(ioHitIndex.filePath());
        }

        //- check and remove basePoints if they are same as current points
        bool samePoints = basepoints == cpoints;
        reduce(samePoints, andOp<bool>());
        if (samePoints)
        {
            rm(basepoints.filePath());
        }
    }
