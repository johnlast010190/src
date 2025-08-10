const dictionary& simple = mesh.solution().dict().subDict("SIMPLE");

const int nNonOrthCorr =
    simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    label tsi =
        simple.lookupOrDefault<label>("thermalSolutionInterval", 1);

    bool therm = true;

    if (tsi > 0)
    {
        therm = (runTime.timeIndex()%tsi == 0);
    }
    else
    {
        therm = false;
    }