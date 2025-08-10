/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow || (!flow && therm)) //flow (U, turb, rho) should be set to steady state
{
    #include "EEqn.H"
}
else if (therm)
{
    if (oCorr == 0)
    {
        if (word(mesh.schemes().ddtScheme("ddt(" + rho.name() + ')')) == "steadyState")
        {
            rho = thermo.rho();
            rho.relax();
        }
        else
        {
            #include "cfdTools/compressible/rhoEqn.H"
        }
    }

    #include "UEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    turb.correct();

    rho = thermo.rho();
    if (word(mesh.schemes().ddtScheme("ddt(" + rho.name() + ')')) == "steadyState")
    {
        rho.relax();
    }
}
else //flow only - unlikely scenario
{
    if (oCorr == 0)
    {
        if (word(mesh.schemes().ddtScheme("ddt(" + rho.name() + ')')) == "steadyState")
        {
            rho = thermo.rho();
            rho.relax();
        }
        else
        {
            #include "cfdTools/compressible/rhoEqn.H"
        }
    }

    #include "UEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    turb.correct();

    rho = thermo.rho();
    if (word(mesh.schemes().ddtScheme("ddt(" + rho.name() + ')')) == "steadyState")
    {
        rho.relax();
    }
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
