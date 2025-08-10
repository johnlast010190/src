
CtrainPtrL_.setSize(phases_.size());
Info<< "create CT fields" << endl;
for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
{
    string phase_Ctname = "Ct";
    string phase_no_str;
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
    phase_Ctname.append(phase_no_str);

    volScalarField Cttemp
    (
        IOobject
        (
            phase_Ctname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(phase_Ctname,dimensionSet(0,0,0,0,0,0,0),1)
    );
    CtrainPtrL_.set(phase_no,Cttemp);
}

if (solveTD_)
{
    if (mesh_.foundObject<volScalarField>("k"))
    {
        // DANGEROUS: can possibly delete the registered object when autoPtr destructor is called
        kPtr_.reset(&const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("k")));
    }
    else
    {
        kPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "k",
                    "0",//runTime.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
        k() *= pow(scalingFactor_,2);
    }
    k().write();
    k().correctBoundaryConditions();

    if (mesh_.foundObject<volScalarField>("epsilon"))
    {
        // DANGEROUS: can possibly delete the registered object when autoPtr destructor is called
        epsilonPtr_.reset(&const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("epsilon")));
    }
    else
    {
        epsilonPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "epsilon",
                    "0",//runTime.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
        epsilon() *= pow(scalingFactor_,3);
    }
    //epsilon().correctBoundaryConditions();
}

if (mesh_.foundObject<volScalarField>("nut"))
{
    // DANGEROUS: can possibly delete the registered object when autoPtr destructor is called
    nutPtr_.reset(&const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("nut")));
}
else
{
    nutPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "nut",
                "0",//runTime.timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
    );
    nut() *= scalingFactor_;
}
nut().correctBoundaryConditions();

nutrainPtr_.reset
(
    new volScalarField
    (
        IOobject
        (
            "nutrain",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("nutrain",dimensionSet(0,2,-1,0,0,0,0),VSMALL)
    )
);
