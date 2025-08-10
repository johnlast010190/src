    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    IOobject meshDict
    (
        dynamicFvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::READ_IF_PRESENT
    );

    bool dynamicMeshExists(false);
    {
        IOobject io
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        dynamicMeshExists = io.typeHeaderOk<IOdictionary>(true);
    }

    autoPtr<dynamicFvMesh> meshPtr;

    if (dynamicMeshExists)
    {
        meshPtr =
            dynamicFvMesh::New
            (
                meshDict
            );
    }
    else
    {
        Info<< "dynamicMeshDict not found" << endl;
        Info<< "Creating staticFvMesh on the fly" << endl;

        meshPtr = new staticFvMesh(meshDict);
    }


    dynamicFvMesh& mesh = meshPtr();
