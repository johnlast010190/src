CorrectPhi
(
    Uc,
    phic,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

// Explicitly choose version of this file supplied by DPMFoam
#include "../DPMFoam/continuityErrs.H"
