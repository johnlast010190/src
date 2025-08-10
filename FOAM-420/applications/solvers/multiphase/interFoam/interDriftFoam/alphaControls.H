const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

scalar blend(alphaControls.lookupOrDefault<scalar>("blend", scalar(0.99)));
