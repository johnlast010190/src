    const bool meshWithChangers =
        runTime.controlDict().lookupOrDefault<bool>("meshChangers", false);

    autoPtr<fvMesh> meshPtr;

    if (meshWithChangers)
    {
        // Create mesh with changers (mesh movers and/or topo changers)
        Info<< "Create mesh with changers for time = "
            << runTime.timeName() << nl << endl;

        meshPtr = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            ),
            true
        );
    }
    else
    {
        Info<< "Create mesh for time = "
            << runTime.timeName() << nl << endl;

        meshPtr = dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        ).ptr();
    }

    fvMesh& mesh = meshPtr();

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
