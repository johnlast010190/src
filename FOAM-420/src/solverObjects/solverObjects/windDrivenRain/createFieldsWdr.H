Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Rh_ = transportProperties.lookup("Rh");
temp_ = transportProperties.lookup("temp");

scalingFactor_ =
(
    transportProperties.lookupOrDefault<scalar>("scalingFactor", 1.0)
);

solveTD_ = transportProperties.lookupOrDefault("solveTD", false);

phases_ = List<List<scalar>>(
    transportProperties.lookup("phases")
);

Info<< "Reading field U\n" << endl;

if (mesh_.foundObject<volVectorField>("U"))
{
    UPtr_ = mesh_.lookupObjectRefPtr<volVectorField>("U");
}
else
{
    UPtr_ =
        new volVectorField
        (
            IOobject
            (
                "U",
                "0",//runTime.timeName(), //If read from the latest timedir, it will be scaled each time the simulation is restarted!
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
    );
    UPtr_->store();
    U() *= scalingFactor_;
    U().correctBoundaryConditions();
}
/*
if (mesh_.foundObject<surfaceScalarField>("phi"))
{
    // DANGEROUS: can possibly delete the registered object when autoPtr destructor is called
    phiwindPtr_.reset(&const_cast<surfaceScalarField&>(mesh_.lookupObject<surfaceScalarField>("phi")));
}
else
{
    phiwindPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiwind",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            linearInterpolate(U()) & mesh.Sf()
        )
    );
}
*/
RePtr_.reset
(
    new volScalarField
    (
        IOobject
        (
            "Re",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Re",dimensionSet(0,0,0,0,0,0,0),500)
    )
);

CdRePtr_.reset
(
    new volScalarField
    (
        IOobject
        (
            "CdRe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("CdRe",dimensionSet(0,0,0,0,0,0,0),250)
    )
);
