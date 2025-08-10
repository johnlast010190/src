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

faMesh aMesh (mesh,faPatch);
