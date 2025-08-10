const dictionary& pimple = mesh.solution().dict().subDict("PIMPLE");

int nNonOrthCorr =
    pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    label tsi =
        pimple.lookupOrDefault<label>("thermalSolutionInterval", 1);

    bool therm = true;

    if (tsi > 0)
    {
        therm = (runTime.timeIndex()%tsi == 0);
    }
    else
    {
        therm = false;
    }
