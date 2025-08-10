Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
thermo.validate(args.executable(), "e");

volScalarField& p = thermo.p();

// Create rho with fixedValueZone boundaries for GIB so that we don't
// try to calculate from continuity in the inactive region
wordList bTypes(p.boundaryField().size(), "calculated");
forAll(p.boundaryField(), patchI)
{
    if (isA<fixedValueZoneFvPatchField<scalar>>(p.boundaryField()[patchI]))
    {
        bTypes[patchI] = fixedValueZoneFvPatchField<scalar>::typeName;
    }
}
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
    thermo.rho(),
    bTypes
);
forAll(rho.boundaryField(), patchI)
{
    if (isA<fixedValueZoneFvPatchField<scalar>>(rho.boundaryField()[patchI]))
    {
        fixedValueZoneFvPatchField<scalar>& fvzpf =
            refCast<fixedValueZoneFvPatchField<scalar>>(rho.boundaryFieldRef()[patchI]);
        fvzpf.cellZoneValue() = gAverage(thermo.rho(patchI));
    }
}


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

mesh.schemes().setFluxRequired(p.name());

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

#include "cfdTools/general/include/createK.H"

#include "cfdTools/general/include/createFvOptions.H"
