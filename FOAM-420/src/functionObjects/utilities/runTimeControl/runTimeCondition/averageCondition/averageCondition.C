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
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeCondition/averageCondition/averageCondition.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(averageCondition, 0);
    addToRunTimeSelectionTable(runTimeCondition, averageCondition, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::averageCondition::averageCondition
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
    fieldNames_(dict.lookup("fields")),
    tolerance_(readScalar(dict.lookup("tolerance"))),
    window_(dict.lookupOrDefault<scalar>("window", -1)),
    startTime_(dict.lookupOrDefault<scalar>("startTime", -1)),
    totalTime_(fieldNames_.size(), obr_.time().deltaTValue()),
    resetOnRestart_(false)
{
    if (resetOnRestart_)
    {
        const dictionary& dict = conditionDict();

        forAll(fieldNames_, fieldi)
        {
            const word& fieldName = fieldNames_[fieldi];

            if (dict.found(fieldName))
            {
                const dictionary& valueDict = dict.subDict(fieldName);
                totalTime_[fieldi] = readScalar(valueDict.lookup("totalTime"));
            }
        }
    }

    //set up file IO
    if (Pstream::master() || !Pstream::parRun())
    {
        //add headers
        writeCommented(file(), "Time");
        forAll(fieldNames_, fieldi)
        {
            word fieldName = fieldNames_[fieldi];
            writeDelimited(file(), fieldName + "-Average");
            writeDelimited(file(), fieldName + "-Delta");
        }
        file() <<endl;
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::averageCondition::~averageCondition()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::averageCondition::apply()
{
    bool satisfied = true;

    if (!active_)
    {
        return satisfied;
    }

    scalar currentTime = obr_.time().value();

    if (currentTime<startTime_)
    {
        return false;
    }

    scalar dt = obr_.time().deltaTValue();

    if (Pstream::master() || !Pstream::parRun())
    {
        file() << obr_.time().timeName();
    }

    Log << "    " << type() << ": " << name_ << " averages:" << nl;

    DynamicList<label> unprocessedFields(fieldNames_.size());

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName(fieldNames_[fieldi]);

        scalar Dt = totalTime_[fieldi];
        scalar alpha = (Dt - dt)/Dt;
        scalar beta = dt/Dt;

        if (window_ > 0)
        {
            if (Dt - dt >= window_)
            {
                alpha = (window_ - dt)/window_;
                beta = dt/window_;
            }
            else
            {
                // Ensure that averaging is performed over window time
                // before condition can be satisfied
                satisfied = false;
            }
        }

        bool processed = false;
        calc<scalar>(fieldName, alpha, beta, satisfied, processed);
        calc<vector>(fieldName, alpha, beta, satisfied, processed);
        calc<sphericalTensor>(fieldName, alpha, beta, satisfied, processed);
        calc<symmTensor>(fieldName, alpha, beta, satisfied, processed);
        calc<tensor>(fieldName, alpha, beta, satisfied, processed);

        if (!processed)
        {
            unprocessedFields.append(fieldi);
        }

        totalTime_[fieldi] += dt;
    }

    if (Pstream::master() || !Pstream::parRun())
    {
        file() << endl;
    }

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        forAll(unprocessedFields, i)
        {
            label fieldi = unprocessedFields[i];
            Info<< "        " << fieldNames_[fieldi] << nl;
        }
    }

    Log << endl;

    return satisfied;
}


void Foam::functionObjects::runTimeControls::averageCondition::write()
{
    dictionary& conditionDict = this->conditionDict();

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        // value dictionary should be present - mean values are written there
        if (conditionDict.found(fieldName))
        {
            dictionary& valueDict = conditionDict.subDict(fieldName);
            valueDict.add("totalTime", totalTime_[fieldi], true);
        }
        else
        {
            dictionary valueDict;
            valueDict.add("totalTime", totalTime_[fieldi], true);
            conditionDict.add(fieldName, valueDict);
        }
    }
}


// ************************************************************************* //
