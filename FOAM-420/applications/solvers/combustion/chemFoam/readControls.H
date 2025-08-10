// #include "cfdTools/general/include/readTimeControls.H"

runTime.controlDict().lookup("adjustTimeStep") >> adjustTimeStep;

maxDeltaT = readScalar(runTime.controlDict().lookup("maxDeltaT"));

minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);