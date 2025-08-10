CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    divU,
    pimple
);

#include "cfdTools/incompressible/continuityErrs.H"
