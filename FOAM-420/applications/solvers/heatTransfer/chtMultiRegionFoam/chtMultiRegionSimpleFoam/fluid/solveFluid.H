/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

//  Pressure-velocity SIMPLE corrector
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);


    if (flow && !frozenFlow)
    {
        p_rgh.storePrevIter();
        rho.storePrevIter();
    }
    {
        if (therm && (frozenFlow || !flow))
        {
            #include "EEqn.H"
        }
        else if (therm && !frozenFlow && flow)
        {
            #include "UEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"

            turb.correct();
        }
        else if (flow)
        {
            #include "UEqn.H"
            #include "pEqn.H"

            turb.correct();
        }
    }
