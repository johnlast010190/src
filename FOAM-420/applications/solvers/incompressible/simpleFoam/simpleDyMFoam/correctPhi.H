CorrectPhi
(
    U,
    phi,
    p,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    geometricZeroField(),
    simple
);

#include "cfdTools/incompressible/continuityErrs.H"
