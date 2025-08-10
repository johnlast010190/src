/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2018, 2020 Esi Ltd
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeCondition/confidenceIntervalCondition/confidenceIntervalCondition.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/etcFiles/etcFiles.H"
#include "primitives/functions/Function1/Table/TableBase.H"
#include "primitives/functions/Function1/Table/Table.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(confidenceIntervalCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        confidenceIntervalCondition,
        dictionary
    );
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypes,
    2
>::names[] =
{
    "relative",
    "absolute"
};

const Foam::NamedEnum
<
    Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypes,
    2
> Foam::functionObjects::runTimeControls
    ::confidenceIntervalCondition::intvTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimeControls::confidenceIntervalCondition::calc
(
    const word& entryName,
    const scalar& deltaT,
    bool& satisfied,
    bool& processed
)
{
    const word valueType =
        state_.objectResultType(functionObjectName_, entryName);

    if (pTraits<scalar>::typeName != valueType)
    {
        return;
    }

    //retrieve index and increment
    const word sampleName(entryName + "NSamplingSteps");
    const word timeName(entryName + "SamplingTime");
    label n = getConditionResult<label>(sampleName, 0);
    scalar t = getConditionResult<scalar>(timeName, 0.0);

    //get current value
    scalar X
        = state_.getObjectResult<scalar>(functionObjectName_, entryName);

    //get mean value
    const word meanName(entryName + "Mean");
    scalar Xmean = getConditionResult<scalar>(meanName, 0.0);

    //calculate difference in X relative to the old mean
    scalar dX = X - Xmean;

    //calculate the new mean and write to state
    if (t != 0.0)
    {
        Xmean = Xmean + deltaT * dX / t;
    }
    setConditionResult(meanName, Xmean);

    //retrieve and update the sum of the square deviation
    const word ssdName(entryName + "SSD");
    scalar XSSD = this->getConditionResult(ssdName, 0.0);
    XSSD += deltaT * dX * (X - Xmean);
    setConditionResult(ssdName, XSSD);

    //update variance
    scalar sigmaSqr = 0;
    if (n > 1)
    {
        sigmaSqr = (scalar(n)/(scalar(n) - 1.0)) * XSSD / t;
    }

    // write standard deviation
    const word stdDevName(entryName + "StdDev");
    setConditionResult(stdDevName, sqrt(sigmaSqr));

    //Update confidence interval
    scalar sigmaErr = sqrt(sigmaSqr/max(1.0, scalar(n)));
    scalar confIntv = Cst_*sigmaErr;
    const word confIntvName(entryName + "Conf");
    setConditionResult(confIntvName, confIntv);

    //calculate tolerance
    scalar tolerance = interval_;

    if (intervalType_ == itRel)
    {
        tolerance *= Xmean;
    }

    Log << entryName << nl
        << tab << "value                       - " << X << nl
        << tab << "N samples                   - " << n << nl
        << tab << "mean                        - " << Xmean << nl
        << tab << "variance                    - " << sigmaSqr << nl
        << tab << "std error                   - " << sigmaErr << nl
        << tab << "confidence level            - "
               << confidence_ * 100.0 << "%" << nl
        << tab << "target confidence interval  - " << tolerance << nl
        << tab << "current confidence interval - " << confIntv << nl;


    if (Pstream::master() || !Pstream::parRun())
    {
        file() << n << tab << t << tab << X << tab << Xmean << tab
               << sigmaSqr << tab << sigmaErr << tab << confIntv
               << tab << Xmean + confIntv << tab << Xmean - confIntv
               << endl;
    }

    setConditionResult(sampleName, n + 1);
    setConditionResult(timeName, t + deltaT);



    //check convergence
    if (confIntv > tolerance || n < minSamples_ || t < minSampleTime_)
    {
        satisfied = false;
    }

    processed = true;
}


Foam::scalar Foam::functionObjects::runTimeControls::confidenceIntervalCondition
::lookupSTD(const scalar& conf) const
{
    IFstream distFile(findEtcFile(tDistFile_));

    dictionary dict(distFile);

    Function1Types::Table<scalar> STD("distribution", dict);

    return STD.value(1-conf);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::confidenceIntervalCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    writeFile(obr_, name, typeName, dict),
    functionObjectName_(dict.lookup("functionObject")),
    entryNames_(dict.lookup("fields")),
    confidence_(readScalar(dict.lookup("confidence"))),
    interval_(readScalar(dict.lookup("targetInterval"))),
    intervalType_
    (
        intvTypeNames_.read
        (
            IStringStream
            (
                dict.lookupOrDefault<word>
                (
                    "intervalType",
                    word("relative")
                )
            )()
        )
    ),
    minSamples_(dict.lookupOrDefault<scalar>("minSamples", 500)),
    minSampleTime_(dict.lookupOrDefault<scalar>("minSampleTime", 0.0)),
    tDistFile_
    (
        dict.lookupOrDefault<fileName>
        (
            "distributionFile",
            "dictData/tDistribution-n10000.dat"
        )
    ),
    Cst_(lookupSTD(confidence_))
{

    // confidence interval has no meaning for fewer than 2 samples
    minSamples_ = max(minSamples_, 2);

    //set up file IO
    if (Pstream::master() || !Pstream::parRun())
    {
        //add headers
        writeCommented(file(), "nSamples");
        writeCommented(file(), "Time");
        writeCommented(file(), "Value");
        writeCommented(file(), "Mean");
        writeCommented(file(), "Variance");
        writeCommented(file(), "Std-error");
        writeCommented(file(), "Range");
        writeCommented(file(), "Upper");
        writeCommented(file(), "Lower");

        file() <<endl;
    }

}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::~confidenceIntervalCondition()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::apply()
{
    bool satisfied = true;

    if (!active_)
    {
        return satisfied;
    }

    scalar dt = obr_.time().deltaTValue();

    DynamicList<label> unprocessedFields(entryNames_.size());

    forAll(entryNames_, entryi)
    {
        const word& entryName(entryNames_[entryi]);

        bool processed = false;

        Log << "    " << type() << ": ";

        calc(entryName, dt, satisfied, processed);

        if (!processed)
        {
            unprocessedFields.append(entryi);
        }
    }

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        forAll(unprocessedFields, i)
        {
            label entryi = unprocessedFields[i];
            Info<< "        " << entryNames_[entryi] << nl;
        }
    }

    Log << endl;

    return satisfied;
}


void Foam::functionObjects::runTimeControls
::confidenceIntervalCondition::write()
{
}


// ************************************************************************* //
