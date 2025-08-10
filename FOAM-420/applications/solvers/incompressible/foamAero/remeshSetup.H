if (firstTimeStep)
{
    U.boundaryFieldRef().updateCoeffs();
    phi = (linearInterpolate(U) & mesh.Sf());
    phi.oldTime();
    U.oldTime();
}

// mesh motion and mesh quality parameters
label nRelaxedErrors = 0;
{
    const dictionary& relaxDict =
    (
        snapDict().found("relaxed")
        ? snapDict().subDict("relaxed")
        : snapDict().subDict("meshQualityControls")
    );

    labelHashSet wrongFaces(mesh.nFaces()/100);
    nRelaxedErrors = motionSmoother::checkMesh(false, mesh, relaxDict, wrongFaces);
}


label nCriticalErrors = 0;
if (snapDict().found("critical"))
{
    const dictionary& criticalDict = snapDict().subDict("critical");

    labelHashSet wrongFaces(mesh.nFaces()/100);
    nCriticalErrors = motionSmoother::checkMesh(false, mesh, criticalDict, wrongFaces);
}

scalar mVel = gMax(mag(U.primitiveField()));

if
(
    (nRelaxedErrors + nCriticalErrors > nErrorsInitial
     && nRelaxedErrors + nCriticalErrors >= allowErrors)
    || nCriticalErrors > nCriticalErrorsInitial
    || (mVel > maxVelocity && mVel > 1.5*initialMaxVelocity)
    || ((runTime.value() - t0) > maxTime)
)
{
    Info<<"Remesh Latest Time: "<<runTime.time().timeOutputValue()<<endl;
    Info<<"Remesh Old Time: "
        <<runTime.time().timeOutputValue() - runTime.time().deltaTValue()
        <<endl;
    Info<<"Remesh Time Step: "<<runTime.time().deltaTValue()<<endl;

    mesh.movePoints(mesh.oldPoints());

    runTime.writeNow();
    runTime.functionObjects().execute();

    label infoI = 0;
    dictionary& geometryDict = snapDict().subDict("geometry");

    forAllConstIter(dictionary, geometryDict, iter)
    {
        const word& key = iter().keyword();
        dictionary& sdict = geometryDict.subDict(key);

        if (sdict.found("transforms"))
        {
            PtrList<dictionary> transforms =
                PtrList<dictionary>(sdict.lookup("transforms"));
            forAll(transforms, tI)
            {
                dictionary& trans = transforms[tI];
                if (trans.found("motionFile"))
                {

                    trans.add("translateVec",transInfo[infoI].first(),true);
                    infoI++;
                }
            }
            sdict.remove("transforms");
            sdict.add("transforms", transforms, false);
        }
    }

    snapDict().regIOobject::writeObject
    (
        runTime.writeFormat(),
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        true
    );

    Pstream::barrier();
    Foam::sleep(5);

    return 1;
}
else
{
    forAll(transInfo, tI)
    {
        vector& tVec = transInfo[tI].first();
        const interpolationTable<vector> iTable = transInfo[tI].second();
        const vector mVel = iTable(runTime.time().timeOutputValue());

        tVec += mVel*runTime.time().deltaTValue();

        if
        (
            mag(tVec.x()) >= maxTranslation.x()
            || mag(tVec.y()) >= maxTranslation.y()
            || mag(tVec.z()) >= maxTranslation.z()
        )
        {
            runTime.writeNow();
            runTime.functionObjects().execute();

            label infoI = 0;
            dictionary& geometryDict = snapDict().subDict("geometry");

            forAllConstIter(dictionary, geometryDict, iter)
            {
                const word& key = iter().keyword();
                dictionary& sdict = geometryDict.subDict(key);

                if (sdict.found("transforms"))
                {
                    PtrList<dictionary> transforms =
                        PtrList<dictionary>(sdict.lookup("transforms"));
                    forAll(transforms, tI)
                    {
                        dictionary& trans = transforms[tI];
                        if (trans.found("motionFile"))
                        {
                            trans.add("translateVec",transInfo[infoI].first(),true);
                            infoI++;
                        }
                    }
                    sdict.remove("transforms");
                    sdict.add("transforms", transforms, false);
                }
            }

            snapDict().regIOobject::writeObject
            (
                runTime.writeFormat(),
                IOstream::currentVersion,
                IOstream::UNCOMPRESSED,
                true
            );

            Pstream::barrier();
            Foam::sleep(5);

            return 1;
        }
    }
}

firstTimeStep =  false;
