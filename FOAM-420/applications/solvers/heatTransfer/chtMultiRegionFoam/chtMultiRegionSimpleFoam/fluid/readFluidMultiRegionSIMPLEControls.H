    const dictionary& simple = fluidRegions[i].solution().dict().subDict("SIMPLE");

    const int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    //solution frequency controls
    label msi =
        simple.lookupOrDefault<label>("flowSolutionInterval", 1);

    label tsi =
        simple.lookupOrDefault<label>("thermalSolutionInterval", 1);

    Switch modifiedMomentumInterp =
        simple.lookupOrDefault<Switch>("modifiedMomentumInterp", false);

    bool flow = true;
    bool therm = true;

    if (msi > 0)
    {
        flow = (runTime.timeIndex()%msi == 0);
    }
    else
    {
        flow = false;
    }

    if (tsi > 0)
    {
        therm = (runTime.timeIndex()%tsi == 0);
    }
    else
    {
        therm = false;
    }
