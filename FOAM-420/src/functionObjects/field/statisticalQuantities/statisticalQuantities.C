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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "statisticalQuantities/statisticalQuantities.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(statisticalQuantities, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        statisticalQuantities,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary& Foam::functionObjects::statisticalQuantities::setDict()
{
    dictionary& propertyDict = this->propertyDict();

    if (!propertyDict.found(name()))
    {
        propertyDict.add(name(), dictionary());
    }

    return propertyDict.subDict(name());
}

void Foam::functionObjects::statisticalQuantities::populateLists()
{
    if
    (
        runFieldScalarLists_.size() == 0
    )
    {
        populate<scalar>(runFieldScalarLists_, scalarMap_);
    }

    if (sDevAver)
    {
        runFieldMeanScalarLists_ = List<DLList<scalar>>(scalarMap_.size());
        forAll(scalarMap_, i)
        {
            const label mapI = scalarMap_[i];
            const word fieldName(fieldNames_[mapI]+ "MovMean");
            dict_.readIfPresent(fieldName, runFieldMeanScalarLists_[i]);
        }
    }

    dict_.readIfPresent("deltaT", runTimeLists_);
}

void Foam::functionObjects::statisticalQuantities::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Value averages");
    writeCommented(os, "Time");
    forAll(fieldNames_, fieldi)
    {
        if (sAver)
        {
            const word meanName(fieldNames_[fieldi] + "MovMean");
            writeDelimited(os, meanName);
        }
        if (sAver2)
        {
            const word mean2Name(fieldNames_[fieldi] + "MovMean2");
            writeDelimited(os, mean2Name);
        }
        if (sDev)
        {
            const word movStdDevName(fieldNames_[fieldi] + "MovStdDev");
            writeDelimited(os, movStdDevName);
        }
        if (sDevAver)
        {
            const word movStdDevAverName(fieldNames_[fieldi] + "MovStdDevAver");
            writeDelimited(os, movStdDevAverName);
        }
        if (sCov)
        {
            const word movCovName(fieldNames_[fieldi] + "MovCoV");
            writeDelimited(os, movCovName);
        }
        if (sCoDav)
        {
            const word movCoDaVName(fieldNames_[fieldi] + "MovCoDaV");
            writeDelimited(os, movCoDaVName);
        }
        if (sCoVav)
        {
            const word movCoVaVName(fieldNames_[fieldi] + "MovCoVaV");
            writeDelimited(os, movCoVaVName);
        }
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::statisticalQuantities::statisticalQuantities
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    dict_(setDict()),
    functionObjectName_("unknown-functionObject"),
    fieldNames_(),
    sAver(dict.lookupOrDefault<Switch>("movAverage", false)),
    sDev(dict.lookupOrDefault<Switch>("movStdDeviation", false)),
    sAver2(dict.lookupOrDefault<Switch>("movAverage2", false)),
    sDevAver(dict.lookupOrDefault<Switch>("movStdDevOfMovAverage", false)),
    sCov(dict.lookupOrDefault<Switch>("CoV", false)),
    sCoDav(dict.lookupOrDefault<Switch>("CoDaV", false)),
    sCoVav(dict.lookupOrDefault<Switch>("CoVaV", false)),
    sampleSize_(dict.lookupOrDefault<label>("sampleSize", 1e+6))
{
    if (sDev)
    {
        sAver = true;
    }
    if (sCoDav)
    {
        sAver = true;
        sAver2 = true;
    }
    if (sCov)
    {
        sAver = true;
        sDev = true;
    }
    if (sCoVav)
    {
        sAver = true;
        sDev = true;
        sDevAver = true;
    }

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::statisticalQuantities::~statisticalQuantities()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::statisticalQuantities::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);
    writeFile::read(dict);

    dict.lookup("functionObject") >> functionObjectName_;
    dict.lookup("fields") >> fieldNames_;

    return true;
}


bool Foam::functionObjects::statisticalQuantities::execute()
{
    if (sampleSize_<=0)
    {
        Log << type() << ": " << name()
            << " sampleSize "
            << sampleSize_
            << " not valid. "
            <<  nl;
        return true;
    }

    //- do it here to leave the other FO to be executed first
    //  Called only once
    populateLists();

    Log << type() << ": " << name()
        << " statistical data using the last existing "
        << sampleSize_ << " data points:" <<  nl;

    file() << time_.timeName();

    calc<scalar>(runFieldScalarLists_, scalarMap_);

    file()<< endl;

    Log << endl;

    return true;
}


bool Foam::functionObjects::statisticalQuantities::write()
{
    return true;
}


// ************************************************************************* //
