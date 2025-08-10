const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
