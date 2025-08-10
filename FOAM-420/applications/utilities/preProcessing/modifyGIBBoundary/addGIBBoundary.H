
    Info<< "Adding GIB boundaries" <<endl;

    // Read patch construct info from dictionary
    PtrList<entry> patchSources(cGIBDict.lookup("boundary"));

    if (patchSources.size())
    {
        // Old and new patches.
        DynamicList<polyPatch*> allPatches(patches.size()+patchSources.size());

        label startFaceI = mesh.nInternalFaces();

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
                const indirectPolyPatch& inpp = refCast<const
                    indirectPolyPatch>(pp);
                allPatches.append
                (
                    inpp.clone
                    (
                        patches
                    ).ptr()
                );
            }
        }
        forAll(patchSources, addedI)
        {
            if (patchSources[addedI].dict().found("faceZoneName"))
            {
                const dictionary& gibDictI = patchSources[addedI].dict();
                gibDictI.lookup("faceZoneName") >> faceZoneName;
            }
            allPatches.append
            (
                polyPatch::New
                (
                    patchSources[addedI].keyword(),
                    patchSources[addedI].dict(),
                    patches.size()+addedI,
                    patches
                ).ptr()
            );
        }
        allPatches.shrink();
        mesh.removeFvBoundary();
        mesh.addFvPatches(allPatches);
    }
    patches.write();
