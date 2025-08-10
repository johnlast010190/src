forAll(e,i)
{
    if (e[i] < eMin)
    {
        e[i] = eMin;
        rhoE[i] = rho[i]*(e[i] + 0.5*magSqr(U[i]));
    }
}

e.correctBoundaryConditions();
