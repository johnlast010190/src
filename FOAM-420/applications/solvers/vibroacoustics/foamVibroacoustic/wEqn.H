wLap = fac::laplacian(w);
w.correctBoundaryConditions();

wLapLap = fac::laplacian(wLap);
//wLapLap.correctBoundaryConditions();

faScalarMatrix wEqn
(
        fam::d2dt2(w)
        + plateD * wLapLap
        + damp1 * fam::ddt(w)
        - damp0 * sqrt(plateD) * fac::ddt(wLap)
        + damp2 * plateD * fac::ddt(wLapLap)
        - load/plateCrossMass
);

wEqn.solve();

w.correctBoundaryConditions();

scalar wMax = gMax(w.primitiveField());
scalar wMin = gMin(w.primitiveField());
Info<< "   wMaxMin = " << wMax << "   " << wMin << endl;
