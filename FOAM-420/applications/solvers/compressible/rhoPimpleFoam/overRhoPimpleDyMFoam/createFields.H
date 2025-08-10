Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "cfdTools/compressible/compressibleCreatePhi.H"

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

mesh.schemes().setFluxRequired(p.name());

#include "cfdTools/compressible/createDpdt.H"

#include "cfdTools/general/include/createK.H"

#include "cfdTools/general/include/createFvOptions.H"


//- Overset specific

// Add solver-specific interpolations
{
    dictionary oversetDict;
    oversetDict.add("U", true);
    oversetDict.add("p", true);
    oversetDict.add("HbyA", true);
    oversetDict.add("grad(p)", true);

    const_cast<dictionary&>
    (
        mesh.schemes().dict()
    ).add
    (
        "oversetInterpolationRequired",
        oversetDict,
        true
    );
}

// Mask field for zeroing out contributions on hole cells
#include "include/createCellMask.H"
