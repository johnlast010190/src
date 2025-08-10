Info<< "\nReading deflection field w" << endl;

areaScalarField w
(
    IOobject
    (
        "w",
        runTime.timeName(),
        mesh,
        //IOobject::MUST_READ,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedScalar("w",dimensionSet(0,1,0,0,0,0,0),0),
    //fixedValueFaPatchScalarField::typeName
    fixedValueZeroGradientFaPatchScalarField::typeName
);

wordList edgeBCs(w.boundaryField().size());

forAll(edgeBCs, eI)
{
    const faPatchScalarField& peI = w.boundaryField()[eI];

    if (isA<fixedValueFaPatchScalarField>(peI))
        edgeBCs[eI] = fixedValueFaPatchScalarField::typeName;
    //else if (isA<fixedValueZeroGradientFaPatchField<scalar>>(peI))
    //    edgeBCs[eI] = zeroGradientFaPatchScalarField::typeName;
    else
        edgeBCs[eI] = zeroGradientFaPatchScalarField::typeName;
}

areaScalarField wLap
(
    IOobject
    (
        "wLap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fac::laplacian(w),
    edgeBCs,
    w.boundaryField().patchTypes()
);

areaScalarField wLapLap
(
    IOobject
    (
        "wLapLap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fac::laplacian(wLap)
);

areaScalarField acc
(
    IOobject
    (
        "acc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fac::d2dt2(w)
);

areaScalarField load
(
    IOobject
    (
        "load",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar("load",dimensionSet(1,-1,-2,0,0,0,0),0)
);

scalarField dload (load.primitiveField().size(), 0.);

if (!planeWave)
{
    fileName acPresFolder;
    if (Pstream::parRun())
    {
        acPresFolder = fileName("..")/".."/"CAA"/"processor"+std::to_string(Pstream::myProcNo());
    }
    else
    {
        acPresFolder  = fileName("..")/"CAA";
    }
    fileName acPresName = "pac_0";

    scalarIOField pAcInp
    (
         IOobject
         (
            acPresName,
            runTime.system(),
            acPresFolder,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        scalarField(0,0.)
    );

    load.primitiveFieldRef() = pAcInp;
}
