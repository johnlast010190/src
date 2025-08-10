Info<< nl << "Calculating specific catch ratio fields\n" << endl;
for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
{
    string phase_scrname = "scr";
    string phase_no_str;
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
    phase_scrname.append(phase_no_str);

    surfaceScalarField normalvel
    (
        IOobject
        (
            "normalvel",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag( (mesh_.Sf()/mesh_.magSf()) & fvc::interpolate(Urain(phase_no)) )
    );

    volScalarField scrtemp
    (
        IOobject
        (
            phase_scrname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::average(normalvel) * alpharain(phase_no)) * ((3600*1E3)/(Rh*phases_[phase_no][1]))
    );
    scrPtrL_.reset(phase_no,scrtemp);
    scr(phase_no).internalField()=scalar(0);
    scr(phase_no).write();
}

Info<< nl << "Calculating catch ratio field\n" << endl;

volScalarField gcr //initialize global catch ratio field
(
    IOobject
    (
        "gcr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    scr(0)*phases_[0][1]
);

for (int phase_no = 1; phase_no < phases_.size(); phase_no++)
{
    gcr += scr(phase_no)*phases_[phase_no][1];
}
gcr.write();
