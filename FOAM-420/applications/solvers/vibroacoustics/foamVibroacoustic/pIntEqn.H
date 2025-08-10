fvScalarMatrix pIntEqn
(
    fvm::d2dt2(pInt)
    - c02 * fvm::laplacian(pInt)
);
pIntEqn.solve();
