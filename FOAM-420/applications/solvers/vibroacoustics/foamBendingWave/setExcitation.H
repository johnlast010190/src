Info<< "Setting excitation" << endl;

scalar curTimeMod = fabs(fmod(runTime.value(),dtLoad)) ;

if (!planeWave)
{
if (curTimeMod < runTime.deltaTValue() || fabs(curTimeMod-dtLoad) < runTime.deltaTValue())
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

    fileName acPresName = "pac_"+std::to_string(itLoad++);

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

    label loadCycles = dtLoad / runTime.deltaTValue();
    dload = (pAcInp - load.primitiveField()) / loadCycles;
}

scalarField& exc = load.primitiveFieldRef();

exc += dload;

}
else
{
scalarField& exc = load.primitiveFieldRef();
exc = sin( 2 * constant::mathematical::pi * freq  *
                                (  runTime.value()
                                 -  1. / c0.value() * ( aMesh.areaCentres() & dimensionedVector("vec" , dimensionSet(0,-1,0,0,0,0,0),vector(1,0,0) ))
                                )
        );
}
