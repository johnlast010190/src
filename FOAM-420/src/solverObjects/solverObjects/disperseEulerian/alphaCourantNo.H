scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

if (mesh.nInternalFaces())
{
    for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
    {
        surfaceScalarField phi(phid(phase_no));

        scalarField sumPhi
        (
           fvc::surfaceSum(mag(phi))().primitiveField()
        );

        alphaCoNum = max(alphaCoNum, 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue());

        meanAlphaCoNum = max(meanAlphaCoNum,
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue());
    }
}

Info<< "Eulerian phases Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;
