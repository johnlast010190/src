if (frameAcceleration.valid())
{
    dimensionedVector fAcc
    (
        "fAcc",
        g.dimensions(),
        frameAcceleration->value(runTime.timeOutputValue())
    );
    gh = volScalarField("gh", ((g - fAcc) & mesh.C()) - ghRef);
    ghf = surfaceScalarField("ghf", ((g - fAcc) & mesh.Cf()) - ghRef);

    // Updating RTPP dict
    runTimeInfoDict().subDict("results").subDict("interFoam").subDict("vector").set
    (
        "frameAcceleration",
        fAcc.value()
    );
}
