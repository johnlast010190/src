
// Limiting function1 types to "constant" and "table"
// The reason is to give to run-time post-processing
// range (min/max). Which is easier from table and constant
// functions.
vector aVectorMax = Zero;
vector aVectorMin = Zero;
scalar aScalarMax = 0;
scalar maxAPlusG = 0;

const word Func1Type(frameAcceleration().type());
if (Func1Type == Function1Types::Table<vector>::typeName)
{
    const Function1Types::Table<vector>& table =
        dynamic_cast<Function1Types::Table<vector>&>(frameAcceleration());
    const vectorField y(table.y());

    const scalarField xComp(y.component(vector::X)());
    const scalarField yComp(y.component(vector::Y)());
    const scalarField zComp(y.component(vector::Z)());
    aVectorMax =
        vector
        (
            y[findMax(xComp)].x(),
            y[findMax(yComp)].y(),
            y[findMax(zComp)].z()
        );
    aVectorMin =
        vector
        (
            y[findMin(xComp)].x(),
            y[findMin(yComp)].y(),
            y[findMin(zComp)].z()
        );
    aScalarMax = mag(y[findMax(mag(y)())]);
    const scalarField magaPlusG(mag(y - g.value()));
    maxAPlusG = magaPlusG[findMax(magaPlusG)];

}
else if (Func1Type != Function1Types::Constant<vector>::typeName)
{
    const vector constValue
    (
        frameAcceleration->value(runTime.timeOutputValue())
    );
    aVectorMax = constValue;
    aVectorMin = constValue;
    aScalarMax = mag(constValue);
    maxAPlusG = mag(constValue - g.value());
}
else
{
    FatalErrorInFunction
        << "frameAcceleration doesn't support: " << nl
        << Func1Type
        << " Function1 type." << nl
        << abort(FatalError);
}

// scalar sub-dict
dictionary dictScalar("scalar");
dictScalar.add("maxFrameAccelerationMagnitude", aScalarMax);
dictScalar.add
(
    "maxFrameAccelerationMagnitudeMinusGravity",
    maxAPlusG
);

// vector sub-dict
dictionary dictVector("vector");
dictVector.add
(
    "frameAcceleration",
    frameAcceleration->value(runTime.timeOutputValue())
);
dictVector.add("g", g.value());
dictVector.add("maxFrameAccelerationComponents", aVectorMax);
dictVector.add("minFrameAccelerationComponents", aVectorMin);

// interFoam sub-dict
dictionary dictInterFoam("interFoam");
dictInterFoam.add("vector", dictVector);
dictInterFoam.add("scalar", dictScalar);

// Results sub-dict
dictionary dictAcceleration("results");
dictAcceleration.add("interFoam", dictInterFoam);

// One more dummy level to store correctly into the IOdictionary
dictionary dictResults("results");
dictResults.add("results", dictAcceleration);

if (runTime.found("runTimeInfo"))
{
    runTimeInfoDict.set
    (
        runTime.lookupObjectRefPtr<IOdictionary>("runTimeInfo")
    );
    if (!runTimeInfoDict().found("results"))
    {
        runTimeInfoDict().set("results", dictionary("results"));
    }
    runTimeInfoDict().subDict("results").set
    (
        "interFoam",
        dictInterFoam
    );
}
else
{
    // Should initialise only when it is not available.
    runTimeInfoDict.set
    (
        new IOdictionary
        (
            IOobject
            (
                "runTimeInfo",
                runTime.timeName(),
                "uniform",
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dictResults
        )
    );
}