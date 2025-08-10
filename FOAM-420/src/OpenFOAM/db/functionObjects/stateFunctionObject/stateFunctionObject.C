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
    (c) 2015-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"
#include "db/Time/Time.H"

const Foam::word Foam::functionObjects::stateFunctionObject::resultsName_ =
    "results";

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::IOdictionary&
Foam::functionObjects::stateFunctionObject::stateDict() const
{
    return time_.functionObjects().stateDict();
}


Foam::IOdictionary& Foam::functionObjects::stateFunctionObject::stateDict()
{
    return const_cast<IOdictionary&>(time_.functionObjects().stateDict());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stateFunctionObject::stateFunctionObject
(
    const word& name,
    const Time& runTime
)
:
    timeFunctionObject(name, runTime)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary& Foam::functionObjects::stateFunctionObject::propertyDict()
{
    IOdictionary& stateDict = this->stateDict();

    if (!stateDict.found(name()))
    {
        stateDict.add(name(), dictionary());
    }

    return stateDict.subDict(name());
}

bool Foam::functionObjects::stateFunctionObject::setTrigger
(
    const label triggeri
)
{
    IOdictionary& stateDict = this->stateDict();

    label oldTriggeri =
        stateDict.lookupOrDefault<label>("triggerIndex", labelMin);

    if (triggeri > oldTriggeri)
    {
        stateDict.set("triggerIndex", triggeri);
        return true;
    }

    return false;
}


Foam::label Foam::functionObjects::stateFunctionObject::getTrigger() const
{
    const IOdictionary& stateDict = this->stateDict();

    return stateDict.lookupOrDefault<label>("triggerIndex", labelMin);
}


bool Foam::functionObjects::stateFunctionObject::foundProperty
(
    const word& entryName
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(name()))
    {
        const dictionary& baseDict = stateDict.subDict(name());
        return baseDict.found(entryName);
    }

    return false;
}


bool Foam::functionObjects::stateFunctionObject::getDict
(
    const word& entryName,
    dictionary& dict
) const
{
    return getObjectDict(name(), entryName, dict);
}


bool Foam::functionObjects::stateFunctionObject::getObjectDict
(
    const word& objectName,
    const word& entryName,
    dictionary& dict
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(objectName))
    {
        const dictionary& baseDict = stateDict.subDict(objectName);
        if (baseDict.found(entryName) && baseDict.isDict(entryName))
        {
            dict = baseDict.subDict(entryName);
            return true;
        }
    }

    return false;
}


Foam::word Foam::functionObjects::stateFunctionObject::resultType
(
    const word& entryName
) const
{
    return objectResultType(name(), entryName);
}


Foam::word Foam::functionObjects::stateFunctionObject::objectResultType
(
    const word& objectName,
    const word& entryName
) const
{
    word result = word::null;
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            forAllConstIter(dictionary, objectDict, iter)
            {
                const dictionary& dict = iter().dict();

                if (dict.found(entryName))
                {
                    return dict.dictName();
                }
            }
        }
    }

    return result;
}


Foam::List<Foam::word>
Foam::functionObjects::stateFunctionObject::objectResultEntries() const
{
    return objectResultEntries(name());
}


Foam::List<Foam::word> Foam::functionObjects::stateFunctionObject::
objectResultEntries
(
    const word& objectName
) const
{
    DynamicList<word> result(2);

    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            forAllConstIter(dictionary, objectDict, iter)
            {
                const dictionary& dict = iter().dict();
                result.append(dict.toc());
            }
        }
    }

    wordList entries;
    entries.transfer(result);

    return entries;
}


// ************************************************************************* //
