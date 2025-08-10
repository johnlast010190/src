forAll(rho,i)
{
    if (rho[i] < rhoMin)
    {
        rho[i] = rhoMin;
    }
}

rho.correctBoundaryConditions();
