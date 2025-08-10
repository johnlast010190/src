    const dictionary& pimple = mesh.solution().dict().subDict("PIMPLE");

    const int nCorr =
        pimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.lookupOrDefault("momentumPredictor", true);

    //solution frequency controls to allow steady state flow operation
    label msi =
        pimple.lookupOrDefault<label>("flowSolutionInterval", 1);

    label tsi =
        pimple.lookupOrDefault<label>("thermalSolutionInterval", 1);

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
