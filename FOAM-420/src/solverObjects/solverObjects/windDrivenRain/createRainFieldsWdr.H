
UrainPtrL_.setSize(phases_.size());
phirainPtrL_.setSize(phases_.size());
alpharainPtrL_.setSize(phases_.size());
scrPtrL_.setSize(phases_.size());

std::stringstream converter;
for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
{
    string phase_Uname = "U";
    string phase_phiname = "phi";
    string phase_alphaname = "alpha";

    string phase_no_str;
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();

    phase_Uname.append(phase_no_str);
    phase_phiname.append(phase_no_str);
    phase_alphaname.append(phase_no_str);

    volVectorField Utemp
    (
        IOobject
        (
            phase_Uname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    UrainPtrL_.set(phase_no,Utemp);

    surfaceScalarField phitemp
    (
        IOobject
        (
            phase_phiname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(Urain(phase_no)) & mesh_.Sf()
    );
    phirainPtrL_.set(phase_no,phitemp);

    volScalarField alphatemp
    (
        IOobject
        (
            phase_alphaname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    alpharainPtrL_.set(phase_no,alphatemp);

}
