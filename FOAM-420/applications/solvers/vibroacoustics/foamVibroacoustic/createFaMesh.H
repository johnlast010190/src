IOdictionary faSol
(
    IOobject
    (
        "faSolution",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word  faPatch (faSol.subDict("VIBROACOUSTIC").lookup("patch"));
label nCycles (faSol.subDict("VIBROACOUSTIC").lookupOrDefault<label>("vibroCycles",1));
scalar  dtLoad  (dimensionedScalar(faSol.subDict("VIBROACOUSTIC").lookup("dtLoad")).value());
label   itLoad(1);
bool planeWave (faSol.subDict("VIBROACOUSTIC").lookupOrDefault<bool>("planeWave",false));

faMesh aMesh (mesh,faPatch);
