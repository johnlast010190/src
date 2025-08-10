{
    volVectorField Udm(alpha.Udm());

    if (alphaHflux)
    {
        Udm += Hvol();
    }

    // compute settling flux direction for interpolation scheme
    tmp<surfaceScalarField> fluxDir;
    fluxDir = fvc::interpolate(Udm) & mesh.Sf();

    // compute face-interpolated velocity (=flux)
    tmp<surfaceVectorField> Uf =
        downwind<vector>(mesh, fluxDir()).interpolate
        (
            Udm
        );
    tmp<surfaceVectorField> Uupw =
        upwind<vector>(mesh, fluxDir()).interpolate
        (
            Udm
        );
    fluxDir.clear();
    tmp<surfaceVectorField> Ulin = linearInterpolate(Udm);
    Uf.ref() *= blend_;
    Uf.ref() += (1.-blend_) * Ulin();

    // blend downwind-upwind interpolation based on Xf
    scalar alpha23 = alpha.alpha23();

    tmp<surfaceScalarField> alphaf = fvc::interpolate(alpha);
    forAll(alphaf(), fI)
    {
        if (alphaf()[fI] < alpha23)
        {
            Uf.ref()[fI] = blend_*Uupw()[fI] + (1.-blend_)*Ulin()[fI];
        }
    }
    forAll(alphaf().boundaryField(), patchI)
    {
        forAll(alphaf().boundaryField()[patchI], fI)
        {
            if (alpha.boundaryField()[patchI][fI] < alpha23)
            {
                Uf.ref().boundaryFieldRef()[patchI][fI] =
                    blend_*Uupw().boundaryField()[patchI][fI]
                  + (1.-blend_)*Ulin().boundaryField()[patchI][fI];
            }
        }
    }
    alphaf.clear();
    Uupw.clear();
    Ulin.clear();

    //Info<< "Computing solid particle flux phiSolid" << endl;
    phir = mesh.Sf() & Uf;

    Uf.clear();
}
