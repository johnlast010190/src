if (adjustTimeStep)
{
    runTime.setDeltaT(min(dtChem, maxDeltaT));
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;

    //exit Foam if deltaT is below user secified minDeltaT
    if (runTime.deltaTValue() < minDeltaT)
    {
        FatalErrorInFunction
            << "minDeltaT = " << minDeltaT << nl
            << "The computed deltaT is lower than the minDeltaT"
            << exit(FatalError);
    }

}
