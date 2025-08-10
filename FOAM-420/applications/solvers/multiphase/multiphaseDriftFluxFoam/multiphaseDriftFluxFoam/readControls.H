PtrDictionary<phase>& phases = mixture.phases();

const dictionary& alphaControls = mesh.solution().solverDict("alpha.*");

bool MULESadvection
(
    alphaControls.lookupOrDefault<Switch>("MULESadvection", false)
);
label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));
bool alphaApplyPrevCorr
(
    alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
);
scalar rSct(alphaControls.lookupOrDefault<scalar>("rSct", scalar(1.0)));
bool alphaDiffusion
(
    alphaControls.lookupOrDefault<Switch>("alphaDiffusion", true)
);
scalar blend_(alphaControls.lookupOrDefault<scalar>("blend", scalar(0.99)));

// MULES Correction
List<tmp<surfaceScalarField>> talphaPhiCorr0(phases.size());

// control addition of drift flux terms
bool alphaHflux(alphaControls.lookupOrDefault<Switch>("alphaHflux", false));
bool momentumDDtHdm
(
    alphaControls.lookupOrDefault<Switch>("momentumDDtHdm", false)
);
bool momentumDivTauDm
(
    alphaControls.lookupOrDefault<Switch>("momentumDivTauDm", false)
);
