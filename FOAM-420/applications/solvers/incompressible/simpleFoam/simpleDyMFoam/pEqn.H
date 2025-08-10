/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
{
    //volScalarField rAU(1.0/UEqn.A());
    rAU = 1.0/UEqn.A();
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    // check difference between U and prevIter-U
    //(U-U.prevIter())().write();

    if (simple.modifiedMomentumInterp())
    {
        fvOptions.makeAbsolute(phi);

        phiHbyA
            += (1 - mesh.solution().equationRelaxationFactor(U.name()))
            *(phi - (fvc::interpolate(U.prevIter())&mesh.Sf()));
    }

    fvOptions.makeRelative(phiHbyA);

    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), fvOptions);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
//Info<<"solving pressure" << endl;
        pEqn.solve();
//            volScalarField ptwo("ptwo", p);
//            ptwo.write();
//Pout<<"finished solving pressure" << endl;
        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }
//Pout<<"min/max(p): " << min(p) << " , " << max(p) << endl;
//Info<<"min/max(p): " << gMin(p) << " , " << gMax(p) << endl;
/*
Info<<"min/max(pEND): " << gMin(p) << " , " << gMax(p) << endl;
forAll(p.internalField(), cellI)
{
    p.internalField()[cellI] += 1e-05;
}
Info<<"now boundaries ..." << endl;
forAll(p.boundaryField(), patchI)
{
    forAll(p.boundaryField()[patchI], faceI)
    {
        p.boundaryField()[patchI][faceI] += 1e-05;
    }
}
Info<<"min/max(pFINAL): " << gMin(p) << " , " << gMax(p) << endl;
*/
    #include "cfdTools/incompressible/continuityErrs.H"
Info<<"relax pressure"<<endl;
    // Explicitly relax pressure for momentum corrector
    p.relax();
//Pout<<"momentum corrector"<<endl;
    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p);
//Pout<<"momentum corrector 0.0.1"<<endl;
//    U -= rAU*fvc::grad(p);
//    U -= rAtU()*fvc::grad(p);
//Pout<<"momentum corrector 0.1"<<endl;
    U.correctBoundaryConditions();
//Info<<"momentum corrector 0.2"<<endl;
    fvOptions.correct(U);
Info<<"re-compute Uf"<<endl;
    // compute Uf based on U and phi
    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }
}
