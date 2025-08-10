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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2012 Esi Ltd.
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "reverseFieldAverage/reverseFieldAverage.H"
#include "fields/volFields/volFields.H"
#include "db/Time/Time.H"
#include "reverseFieldAverage/reverseFieldAverageItem/reverseFieldAverageItem.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reverseFieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, reverseFieldAverage, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reverseFieldAverage::resetFields()
{
    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            if (obr().found(faItems_[i].meanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[i].meanFieldName()]);
            }
        }

        if (faItems_[i].prime2Mean())
        {
            if (obr().found(faItems_[i].prime2MeanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[i].prime2MeanFieldName()]);
            }
        }
    }
}


void Foam::functionObjects::reverseFieldAverage::initialize()
{
    resetFields();

    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr_.time().deltaTValue());

    // Add mean fields to the field lists
    forAll(faItems_, fieldI)
    {
        addMeanField<scalar>(fieldI);
        addMeanField<vector>(fieldI);
        addMeanField<sphericalTensor>(fieldI);
        addMeanField<symmTensor>(fieldI);
        addMeanField<tensor>(fieldI);
    }

    // Add prime-squared mean fields to the field lists
    forAll(faItems_, fieldI)
    {
        addPrime2MeanField<scalar, scalar>(fieldI);
        addPrime2MeanField<vector, symmTensor>(fieldI);
    }

    forAll(faItems_, fieldI)
    {
        if (!faItems_[fieldI].active())
        {
            WarningInFunction
                << "Field " << faItems_[fieldI].fieldName()
                << " not found in database for averaging";
        }
    }

    // ensure first averaging works unconditionally
    prevTimeIndex_ = -1;

    Log << endl;

    initialised_ = true;
}


void Foam::functionObjects::reverseFieldAverage::restart()
{
    Log  << "    Restarting averaging at time " << obr_.time().timeName()
         << nl << endl;

    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr_.time().deltaTValue());

    initialize();
}


void Foam::functionObjects::reverseFieldAverage::calcAverages()
{
    if (!initialised_)
    {
        initialize();
    }

    const label currentTimeIndex = obr_.time().timeIndex();
    const scalar currentTime = obr_.time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }

    if (periodicRestart_ && currentTime > restartPeriod_*periodIndex_)
    {
        restart();
        periodIndex_++;
    }

    Log << type() << " " << name() << " execute:" << nl
        << "    Calculating averages" << nl;

    addMeanSqrToPrime2Mean<scalar, scalar>();
    addMeanSqrToPrime2Mean<vector, symmTensor>();

    calculateMeanFields<scalar>();
    calculateMeanFields<vector>();
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculatePrime2MeanFields<scalar, scalar>();
    calculatePrime2MeanFields<vector, symmTensor>();

    forAll(faItems_, fieldI)
    {
        totalIter_[fieldI]++;
        totalTime_[fieldI] += obr_.time().deltaTValue();
    }
}


void Foam::functionObjects::reverseFieldAverage::writeAverages() const
{
    Log << "    Writing average fields" << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();
}


void Foam::functionObjects::reverseFieldAverage::writeAveragingProperties()
{
    forAll(faItems_, fieldI)
    {
        const word& fieldName = faItems_[fieldI].fieldName();

        dictionary propsDict;
        propsDict.add("totalIter", totalIter_[fieldI]);
        propsDict.add("totalTime", totalTime_[fieldI]);
        setProperty(fieldName, propsDict);
    }
}


void Foam::functionObjects::reverseFieldAverage::readAveragingProperties()
{
    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr_.time().deltaTValue());

    if (restartOnRestart_ || restartOnOutput_)
    {
        Info<< "    Starting averaging at time " << obr_.time().timeName()
            << nl;
    }
    else
    {
        Info<< "    Restarting averaging for fields:" << nl;

        forAll(faItems_, fieldI)
        {
            const word& fieldName = faItems_[fieldI].fieldName();
            if (foundProperty(fieldName))
            {
                dictionary fieldDict;
                getProperty(fieldName, fieldDict);

                totalIter_[fieldI] = readLabel(fieldDict.lookup("totalIter"));
                totalTime_[fieldI] = readScalar(fieldDict.lookup("totalTime"));

                Info<< "        " << fieldName
                    << " iters = " << totalIter_[fieldI]
                    << " time = " << totalTime_[fieldI] << nl;
            }
            else
            {
                Info<< "        " << fieldName
                    << ": starting averaging at time "
                    << obr_.time().timeName() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reverseFieldAverage::reverseFieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(GREAT),
    averagingStartTime_(runTime.endTime().value()),
    initialised_(false),
    faItems_(),
    totalIter_(0),
    totalTime_(0),
    periodIndex_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::reverseFieldAverage::~reverseFieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reverseFieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    initialised_ = false;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput", restartOnOutput_);
    dict.readIfPresent("periodicRestart", periodicRestart_);
    dict.readIfPresent("averagingStartTime", averagingStartTime_);

    dict.lookup("fields") >> faItems_;
    forAll(faItems_, faI)
    {
        faItems_[faI].setMeanFieldName
        (
            faItems_[faI].meanFieldName() + "_" + name()
        );
        faItems_[faI].setPrime2MeanFieldName
        (
            faItems_[faI].prime2MeanFieldName() + "_" + name()
        );
    }

    if (periodicRestart_)
    {
        dict.lookup("restartPeriod") >> restartPeriod_;
    }

    scalar currentTime = obr_.time().value();
    scalar halfDeltaT = 0.5*obr_.time().deltaTValue();

    if (currentTime > averagingStartTime_ + halfDeltaT)
    {
        readAveragingProperties();
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::reverseFieldAverage::execute()
{
    scalar currentTime = obr_.time().value();
    scalar halfDeltaT = 0.5*obr_.time().deltaTValue();

    if
    (
        currentTime < (averagingStartTime_ + halfDeltaT)
     && prevTimeIndex_ > obr_.time().timeIndex()
    )
    {
        calcAverages();
        if (obr_.time().outputTime())
            write();

        Log << endl;
    }
    prevTimeIndex_ = obr_.time().timeIndex();

    return true;
}


bool Foam::functionObjects::reverseFieldAverage::write()
{
    scalar cTime = obr_.time().value();
    scalar hdeltaT = 0.5*obr_.time().deltaTValue();

    if
    (
        cTime <= (averagingStartTime_ + hdeltaT)
     && prevTimeIndex_ > obr_.time().timeIndex()
    )
    {
        writeAverages();
        writeAveragingProperties();

        if (restartOnOutput_)
        {
            restart();
        }

        Log << endl;
    }

    return true;
}

// ************************************************************************* //
