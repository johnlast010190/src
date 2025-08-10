{
    wordList pcorrTypes
    (
        U.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", dimensionSet(1,-1,-2,0,0), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/dimensionSet(1,-3,0,0,0),
        runTime.deltaT().value()
    );

    adjustPhi(phi, U, pcorr);

    mesh.schemes().setFluxRequired(pcorr.name());

    label nNonOrthCorrs = 3;
    for (label count=1; count <= nNonOrthCorrs; count++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(0, 0);
        pcorrEqn.solve();

        if (count == nNonOrthCorrs)
        {
            phi -= pcorrEqn.flux();
        }
    }
}
