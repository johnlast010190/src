//whether to remesh when quality issues develop
Switch reMesh
(
    pimple.dict().lookupOrDefault<Switch>("foamRemesh",false)
);

vector maxTranslation
(
    pimple.dict().lookupOrDefault<vector>("maxTranslation",vector(GREAT, GREAT, GREAT))
);

label allowErrors
(
    pimple.dict().lookupOrDefault<label>("allowErrors",0)
);

scalar maxVelocity
(
    pimple.dict().lookupOrDefault<scalar>("maxVelocity",GREAT)
);

scalar maxTime
(
    pimple.dict().lookupOrDefault<scalar>("maxTime",GREAT)
);

scalar scaleFirstTimeStep
(
    pimple.dict().lookupOrDefault<scalar>("scaleFirstTimeStep",1.0)
);

autoPtr<IOdictionary> snapDict;
label nErrorsInitial = 0;
label nCriticalErrorsInitial = 0;
scalar initialMaxVelocity = gMax(mag(U.primitiveField()));
scalar t0 = runTime.value();

DynamicList<Tuple2<vector, interpolationTable<vector>>>
transInfo(mesh.boundary().size());

bool firstTimeStep = true;

if (reMesh)
{
    // Read meshing dictionary
    snapDict.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "foamHexMeshDict",
                runTime.caseSystem(),
                runTime,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
             )
         )
     );

    // mesh motion and mesh quality parameters
    {
        const dictionary& motionDict =
            snapDict().subDict("meshQualityControls");

        labelHashSet wrongFacesInit(mesh.nFaces()/100);
        nErrorsInitial = motionSmoother::checkMesh(false, mesh, motionDict, wrongFacesInit);
    }

    if (snapDict().found("critical"))
    {
        const dictionary& criticalDict = snapDict().subDict("critical");

        labelHashSet wrongFaces(mesh.nFaces()/100);
        nCriticalErrorsInitial = motionSmoother::checkMesh(false, mesh, criticalDict, wrongFaces);
    }

    dictionary& geometryDict = snapDict().subDict("geometry");

    DynamicList<Tuple2<vector, interpolationTable<vector>>>
        tInfo(geometryDict.size());

    forAllConstIter(dictionary, geometryDict, iter)
    {
        const word& key = iter().keyword();
        const dictionary& sdict = geometryDict.subDict(key);

        if (sdict.found("transforms"))
        {
            PtrList<dictionary> transforms(sdict.lookup("transforms"));
                //
            forAll(transforms, tI)
            {
                const dictionary& transformDict = transforms[tI];
                if (transformDict.found("motionFile"))
                {
                    vector transVec = transformDict.lookup("translateVec");
                    fileName motionFile = transformDict.lookup("motionFile");
                    interpolationTable<vector> iTable(motionFile);

                    transInfo.append
                    (
                        Tuple2<vector, interpolationTable<vector>>
                        (transVec, iTable)
                     );
                    Info<<"Surface: "<<key<<" Initial Translation vector: "<<transVec<<endl;
                }
            }
        }
    }
}
transInfo.shrink();
