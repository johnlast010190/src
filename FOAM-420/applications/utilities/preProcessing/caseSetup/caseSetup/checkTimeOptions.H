//
// checkTimeOptions.H
// ~~~~~~~~~~~~~~~~~~

instantList Times = runTime.times();

Foam::scalar startTime = 0;

// check -time and -latestTime options
if (args.optionFound("time"))
{
    startTime = args.optionRead<scalar>("time");
}

if (args.optionFound("latestTime"))
{
    if (Times[Times.size() - 1].name() != "constant")
    {
        startTime = Times[Times.size() - 1].value();
    }
}

runTime.setTime
(
    dimensionedScalar("time", dimTime, startTime),
    0
);
