
Switch solveThermal =
    pimple.dict().lookupOrDefault<Switch>("solveThermal", true);

Switch ddtPhiCorr =
    pimple.dict().lookupOrDefault<Switch>("ddtPhiCorr", true);
