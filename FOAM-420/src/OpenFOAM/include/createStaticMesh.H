const bool meshWithChangers =
    runTime.controlDict().lookupOrDefault<bool>("meshChangers", false);

if (meshWithChangers)
{
    // Create static mesh with dummy changers. Make the solver compatible with
    // static non-conformal couplings.
    Foam::Info<< "Create static mesh with changers for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}
else
{
    Foam::Info<< "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    meshWithChangers
);

// Test for the presence of non-conformal patches in a mesh without changers
if (!meshWithChangers)
{
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if (isA<Foam::nonConformalPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Non-conformal patches are not compatible with old "
                << "dynamic mesh structure, without changers." << nl
                << "Please set the 'meshChangers' flag to 'true' "
                << "in the controlDict dictionary or remove any "
                << "non-conformal patch."
                << exit(FatalError);
        }
    }
}
