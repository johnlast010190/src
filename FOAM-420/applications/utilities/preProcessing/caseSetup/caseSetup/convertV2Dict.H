/****************************
Reads caseSetupDict and checks
if this is the old or new version.
if the old version, this will
convert and output the new version
*****************************/



if (!dict.found("global"))
{
    Info<< "Found caseSetupDict version 2: converting to version 3"
         << nl << endl;
    Info<< "Existing file will be saved as <inpufileName>.v2" << endl;

    //save old caseSetup dictionary
    //if (distributed || !Pstream::parRun() || Pstream::master())
    {
        cp(dict.filePath(), dict.filePath() + ".v2");
    }

    const dictionary& v2r0(dict.subDict("region0"));

    // generate v3 setup file
    dictionary v3SetupDict;

    // global
    {
        v3SetupDict.add("global", dictionary(), false);
        dictionary& global = v3SetupDict.subDict("global");

        //add modification switches to global
        if (v2r0.found("modificationSwitches"))
        {
            global.add
            (
                "modificationSwitches",
                v2r0.subDict("modificationSwitches")
            );
        }
        else
        {
            global.add
            (
                "modificationSwitches",
                dictionary()
            );
        }

        if (global.subDict("modificationSwitches").found("reuseExistingDicts"))
        {
            //reuseExistingDicts no longer supported
            Warning << "Reuse of existing dictionaries via the "
                    << "'reuseExistingDicts' switch is no longer supported. "
                    << endl;
            global.subDict("modificationSwitches").remove("reuseExistingDicts");
        }

        //add default region groups
        HashTable<wordList> regionGroups;
        regionGroups.insert
        (
            "fluid",
            wordList(1, polyMesh::defaultRegion)
        );
        global.add("regionGroups", regionGroups);

        //add the system subdict
        global.add("system", dictionary(), false);

        //add the function subdict
        global.add("functions", dictionary(), false);

        //if function list found, add to dict
        if (v2r0.found("functions"))
        {
            PtrList<entry> functionDicts
            (
                v2r0.lookup("functions")
            );

            forAll(functionDicts, i)
            {
                global.subDict("functions").add
                (
                    functionDicts[i]
                );
            }
        }
    }




    //base regions dict
    v3SetupDict.add("regions", dictionary());

    //add fluid group
    v3SetupDict.subDict("regions").add("fluid", dictionary());

    // region0
    v3SetupDict.subDict("regions").add("region0", dictionary());
    dictionary& reg0(v3SetupDict.subDict("regions").subDict("region0"));

    //modification switches
    reg0.add("modificationSwitches", dictionary());

    //boundary mesh
    //we might need to load the mesh, create storage for it
    //autoPtr<fvMesh> tmpMeshPtr;
    {
        //convert partialNamed to explicit
        reg0.add
        (
            "boundaryMesh", dictionary()
        );

        dictionary& bmesh(reg0.subDict("boundaryMesh"));

        if (v2r0.found("boundaryMesh"))
        {

            if (v2r0.subDict("boundaryMesh").found("partialNamed"))
            {
                //load mesh if necessary
                //createMesh(runTime, tmpMeshPtr);
                const fvPatchList& patches = tmpMeshPtr().boundary();

                const dictionary& pnds
                (
                    v2r0.subDict("boundaryMesh").subDict("partialNamed")
                );

                forAllConstIter(dictionary, pnds, iter)
                {
                    word substring(iter().keyword());
                    const dictionary& ed(iter().dict());

                    for (label patchI = 0; patchI < patches.size(); patchI++)
                    {
                        const word& pname = patches[patchI].name();

                        if (pname.find(substring, 0) != string::npos)
                        {
                            bmesh.add(pname, ed);
                        }
                    }
                }
            }

            if (v2r0.subDict("boundaryMesh").found("exactNamed"))
            {
                bmesh.merge(v2r0.subDict("boundaryMesh").subDict("exactNamed"));
            }
        }

    }

    //merge old global into region0
    reg0.merge(v2r0.subDict("global"));

    //add global
    reg0.add("global", dictionary());

    //check for fvOptions and change to new format if found
    if (reg0.found("system"))
    {
        if (reg0.subDict("system").found("fvOptions"))
        {
            dictionary& fvOptions(reg0.subDict("system").subDict("fvOptions"));

            forAllIter(dictionary, fvOptions, iter)
            {
                if (iter().isDict())
                {
                    dictionary& cOpt(iter().dict());

                    word fvoType(cOpt.lookup("type"));
                    word fvoTypeCoeffs(fvoType+"Coeffs");

                    if (cOpt.found(fvoTypeCoeffs))
                    {
                        dictionary& coeffDict(cOpt.subDict(fvoTypeCoeffs));

                        //grab selectionMode and data input
                        word entryName = "selectionMode";

                        if (cOpt.found(entryName))
                        {
                            word sMode(cOpt.lookup(entryName));
                            updateFvOptionEntry(cOpt, coeffDict, sMode);
                        }

                        updateFvOptionEntry(cOpt, coeffDict, entryName);
                    }
                }
            }
        }
    }


    //fields
    {
        //merge fields
        if (v2r0.found("fields"))
        {
            reg0.add("fields", v2r0.subDict("fields"));
        }


        // add boundary type defaults
        const dictionary& v2bcs(v2r0.subDict("boundaryConditions"));
        if (v2bcs.found("regionDefaults"))
        {
            reg0.add
            (
                "boundaryTypeDefaults",
                v2bcs.subDict("regionDefaults")
            );
        }
        else
        {
            reg0.add
            (
                "boundaryTypeDefaults",
                dictionary()
            );
        }


        // Add boundary conditions (groups, wordRe and explicit)
        reg0.add("boundaryConditions", dictionary());
        dictionary& reg0bcs(reg0.subDict("boundaryConditions"));

        // convert partial named to explicit
        if (v2bcs.found("partialNamed"))
        {
            if (v2bcs.subDict("partialNamed").size())
            {
                //createMesh(runTime, tmpMeshPtr);
                const fvPatchList& patches = tmpMeshPtr().boundary();

                forAllConstIter(dictionary, v2bcs.subDict("partialNamed"), iter)
                {
                    word substring(iter().keyword());
                    const dictionary& ed(iter().dict());

                    for (label patchI = 0; patchI < patches.size(); patchI++)
                    {
                        const word& pname = patches[patchI].name();

                        if (pname.find(substring, 0) != string::npos)
                        {
                            reg0bcs.add(pname, ed);
                        }
                    }
                }

            }
        }

        if (v2bcs.found("exactNamed"))
        {
            reg0bcs.merge(v2bcs.subDict("exactNamed"));
        }
    }

    //add empty functions
    reg0.add("functions", dictionary());

    dict.dictionary::operator=(v3SetupDict);

    if (distributed || !Pstream::parRun() || Pstream::master())
    {
        Info<< "Writing version 3 format caseSetup input file to: "
             << dict.filePath() << "... ";

        dict.regIOobject::writeObject
        (
            runTime.writeFormat(),
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
        );
        Info<< "done"  << nl << endl;
    }

}
