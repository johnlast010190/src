/*---------------------------------------------------------------------------*\
| Modified 2010-2012 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

    autoPtr<IOdictionary> dictPtr;
    word dictName("caseSetupDict");

    IOobject caseSetupHeader
    (
        dictName,
        runTime.caseSystem(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (args.optionFound("dict"))
    {
        fileName customDict = args["dict"];

        caseSetupHeader = IOobject
        (
            customDict,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (Pstream::parRun())
        {
            caseSetupHeader.instance() = runTime.caseSystem();
        }

        //Info<< customDict << endl;

        //check for parallel operations


        //Info<< caseSetupHeader.caseName() << endl;
        //Info<< caseSetupHeader.path() << endl;
        //Info<< caseSetupHeader.rootPath() << endl;
        //Info<< caseSetupHeader.instance() << endl;
        //Info<< caseSetupHeader.local() << endl;

    }

    if (!caseSetupHeader.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorIn(args.executable())
            << "Cannot open caseSetup input file\n    "
            << caseSetupHeader.objectPath()
            << nl
            << exit(FatalError);
    }


    //v2 caseSetupDict conversion, only on master
    {
        //check if mesh will be required for conversion and load it if yes
        bool meshRequiredForUpdate = false;

        if
        (
            collated || distributed || Pstream::master() || !Pstream::parRun()
        )
        {
            IOdictionary::globalIOdictionary = 0;
            entry::disableFunctionEntries = 1;

            IOdictionary dict(caseSetupHeader);

            if (!dict.found("global"))
            {
                const dictionary& v2r0(dict.subDict("region0"));
                if (v2r0.found("boundaryMesh"))
                {
                    if (v2r0.subDict("boundaryMesh").found("partialNamed"))
                    {
                        meshRequiredForUpdate = true;
                    }
                }
                const dictionary& v2bcs(v2r0.subDict("boundaryConditions"));
                if (v2bcs.found("partialNamed"))
                {
                    if (v2bcs.subDict("partialNamed").size())
                    {
                        meshRequiredForUpdate = true;
                    }
                }
            }

            entry::disableFunctionEntries = 0;
            IOdictionary::globalIOdictionary = 1;
        }

        autoPtr<fvMesh> tmpMeshPtr;
        reduce(meshRequiredForUpdate, orOp<bool>());
        if (meshRequiredForUpdate)
        {
            createMesh(runTime, tmpMeshPtr);
        }

        if (collated || distributed || Pstream::master() || !Pstream::parRun())
        {
            IOdictionary::globalIOdictionary = 0;
            entry::disableFunctionEntries = 1;

            IOdictionary dict(caseSetupHeader);

    #include "convertV2Dict.H" //convert the old format to the new format
            entry::disableFunctionEntries = 0;
            IOdictionary::globalIOdictionary = 1;
        }
    }



    IOdictionary dict(caseSetupHeader);

    // find defaults
    dictionary defaults;

    if (dict.found("defaults")) //defaults are included in caseSetupDict
    {
        defaults = dict.subDict("defaults");
        dict.remove("defaults");
    }
    else //defaults sourced from caseDicts
    {
        fileNameList defaultDirs
        (
            findEtcFiles("caseDicts/preProcessing/caseSetup/caseSetup.cfg", true)
        );

        forAll(defaultDirs, di)
        {
            defaultDirs[di].expand();
            defaultDirs[di].toAbsolute();

            defaults.merge
            (
                IOdictionary
                (
                    IOobject
                    (
                        defaultDirs[di],
                        runTime,
                        IOobject::MUST_READ
                    )
                ).subDict("defaults")
            );
        }
    }

    // Now that defaults are extracted, replace dict with the un-evaluated
    // version if requested
    if
    (
        dict.subOrEmptyDict("global").subOrEmptyDict("modificationSwitches")
       .lookupOrDefault("preserveFunctionEntries", false)
    )
    {
        int oldFlag = entry::disableFunctionEntries;
        entry::disableFunctionEntries = 1;
        // Need to use local dictionary here as UIPstream does not read the
        // function entries correctly when distributing
        localIOdictionary dict0(caseSetupHeader);
        if (dict0.found("defaults"))
        {
            dict0.remove("defaults");
        }
        // Prevent re-reading file due to implicit creation of IOdictionary
        dict = dict0;
        entry::disableFunctionEntries = oldFlag;
    }
