/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
//    turbulence->R()
    turbulence->devReff()
);

runTime.write();

if (graphData)
{
    const word& gFormat = runTime.graphFormat();

    makeGraph(y, flowDirection & U, "Uf", gFormat);

    makeGraph(y, turbulence->nu(), gFormat);
    makeGraph(y, turbulence->nut(), gFormat);
    makeGraph(y, turbulence->k(), gFormat);
    makeGraph(y, turbulence->epsilon(), gFormat);

    makeGraph(y, flowDirection & R & flowDirection, "Rff", gFormat);
    makeGraph(y, wallNormal12 & R & wallNormal12, "Rww", gFormat);
    makeGraph(y, flowDirection & R & wallNormal12, "Rfw", gFormat);

    makeGraph(y, sqrt(mag(R.component(symmTensor::XX))), "u", gFormat);
    makeGraph(y, sqrt(mag(R.component(symmTensor::YY))), "v", gFormat);
    makeGraph(y, sqrt(mag(R.component(symmTensor::ZZ))), "w", gFormat);
    makeGraph(y, R.component(symmTensor::XY), "uv", gFormat);

    makeGraph(y, mag(fvc::grad(U)), "gammaDot", gFormat);

    scalar tauw = (flowDirection
         & R.boundaryField()[patchId][faceId]
         & wallNormal12);

    dimensionedScalar uTau
    (
        "uTauw",
        dimVelocity,
        ::sqrt(mag(tauw))
    );

    dimensionedScalar nuSMALL
    (
        "nuSMALL",
        dimensionSet(0, 2, -1, 0, 0),
        ROOTVSMALL
    );
    volScalarField tauF
    (
        IOobject
        (
            "tauF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, tauw)
    );
    makeGraph
    (
        y*uTau/(turbulence->nu() + nuSMALL),
        tauF,
        "t_w",
        gFormat
    );

    makeGraph
    (
        y*uTau/(turbulence->nu() + nuSMALL),
        mag(U)/uTau,
        "u+",
        gFormat
    );

    makeGraph
    (
        y*uTau/(turbulence->nu() + nuSMALL),
        turbulence->nut(),
        "nut+",
        gFormat
    );
    makeGraph
    (
        y,
        turbulence->nut(),
        "nut",
        gFormat
    );
    volScalarField omega(turbulence->epsilon()/turbulence->k()/0.09);
    makeGraph
    (
        y*uTau/(turbulence->nu()() + nuSMALL),
        -fvc::grad(omega)&wallNormal12,
        "gradOmega+",
        gFormat
    );

    makeGraph
    (
        y*uTau/(turbulence->nu()() + nuSMALL),
        -fvc::grad(turbulence->k())&wallNormal12,
        "gradK+",
        gFormat
    );
    makeGraph
    (
        y*uTau/(turbulence->nu()() + nuSMALL),
        turbulence->epsilon()/turbulence->k()/0.09,
        "omega_y+",
        gFormat
    );
    makeGraph
    (
        y,
        turbulence->epsilon()/turbulence->k()/0.09,
        "omega_y",
        gFormat
    );
    makeGraph
    (
        y*uTau/(turbulence->nu()() + nuSMALL),
        turbulence->k(),
        "k_y+",
        gFormat
    );
}
