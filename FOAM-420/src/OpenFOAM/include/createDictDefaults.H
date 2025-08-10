/*---------------------------------------------------------------------------*\
| Copyright (C) 2016 Esi Ltd                                               |
\*---------------------------------------------------------------------------*/

Info<< endl;

Switch distributed = false;

if (Pstream::parRun() && Foam::fileHandler().type() == "uncollated")
{
    fileName dpdName("system/decomposeParDict");

    IFstream dds(dpdName);

    if (dds)
    {
        dictionary decomposeParDict(dds);

        if  (decomposeParDict.found("distributed"))
        {
           distributed = readBool(decomposeParDict.lookup("distributed"));
        }
    }

    if (distributed)
    {
        Info<< "Distributed mode: on" << endl;
    }
}

//check if required dicitonaries: controlDict, fvSchemes and fvSolution exist
//if not copy from defaults
if (distributed || !Pstream::parRun() || Pstream::master())
{
    fileName ctrlName("system/controlDict");

    IFstream cs(ctrlName);
    if (!cs)
    {
        Info<< "Creating 'system/controlDict'" << endl;

        cp
        (
            findEtcFile
            (
                "caseDicts/preProcessing/caseSetup/settings/"
                "defaultControlDict.cfg"
            ),
            "system/controlDict"
        );
    }

    word regionName;
    fileName schName;
    fileName solName;

    if (args.optionReadIfPresent("region", regionName))
    {
        schName = "system/" + regionName + "/fvSchemes";
        solName = "system/" + regionName + "/fvSolution";
    }
    else
    {
        schName = "system/fvSchemes";
        solName = "system/fvSolution";
    }

    IFstream scs(schName);
    if (!scs)
    {
        Info<< "Creating " << schName <<endl;

        cp
        (
            findEtcFile
            (
                "caseDicts/preProcessing/caseSetup/settings/"
                "defaultFvSchemes.cfg"
            ),
            schName
        );
    }

    IFstream sos(solName);
    if (!sos)
    {
        Info<< "Creating " << solName <<endl;

        cp
        (
            findEtcFile
            (
                "caseDicts/preProcessing/caseSetup/settings/"
                "defaultFvSolution.cfg"
            ),
            solName
        );
    }
}
