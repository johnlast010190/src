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

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "global/debug/debug.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "global/etcFiles/etcFiles.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "include/demandDrivenData.H"
#include "global/debug/simpleObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! \cond ignoreDocumentation
//- Skip documentation : local scope only

dictionary* controlDictPtr_(nullptr);
dictionary* debugSwitchesPtr_(nullptr);
dictionary* infoSwitchesPtr_(nullptr);
dictionary* optimisationSwitchesPtr_(nullptr);
dictionary* tolerancesPtr_(nullptr);

// Debug switch read and write callback tables.
simpleObjectRegistry* debugObjectsPtr_(nullptr);
simpleObjectRegistry* infoObjectsPtr_(nullptr);
simpleObjectRegistry* optimisationObjectsPtr_(nullptr);
simpleObjectRegistry* dimensionSetObjectsPtr_(nullptr);
simpleObjectRegistry* dimensionedConstantObjectsPtr_(nullptr);


// To ensure controlDictPtr_ is deleted at the end of the run
class deleteControlDictPtr
{
public:

    deleteControlDictPtr()
    {}

    ~deleteControlDictPtr()
    {
        deleteDemandDrivenData(debugObjectsPtr_);
        deleteDemandDrivenData(infoObjectsPtr_);
        deleteDemandDrivenData(optimisationObjectsPtr_);
        deleteDemandDrivenData(dimensionSetObjectsPtr_);
        deleteDemandDrivenData(dimensionedConstantObjectsPtr_);

        debugSwitchesPtr_ = nullptr;
        infoSwitchesPtr_ = nullptr;
        optimisationSwitchesPtr_ = nullptr;
        deleteDemandDrivenData(controlDictPtr_);
    }
};

deleteControlDictPtr deleteControlDictPtr_;
//! \endcond


} // End namespace debug
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        controlDictPtr_ = new dictionary();
        forAllReverse(controlDictFiles, cdfi)
        {
            IFstream ifs(controlDictFiles[cdfi]);

            if (!ifs.good())
            {
                SafeFatalIOErrorInFunction
                (
                    ifs,
                    "Cannot open controlDict"
                );
            }
            controlDictPtr_->merge(dictionary(ifs));
        }
    }

    return *controlDictPtr_;
}


/*
//Non IFstream construction for static binary
Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        dictionary dict = dictionary();

        dict.add(word("Documentation"), dictionary(), false);
        dict.add(word("DebugSwitches"), dictionary(), false);
        dict.add(word("InfoSwitches"), dictionary(), false);
        dict.add(word("OptimisationSwitches"), dictionary(), false);
        char defaultDiv[]
            = "nonBlocking";
        dict.subDict("OptimisationSwitches").add
        (
            word("commsType"),
            defaultDiv,
            true
         );

        dict.add(word("DimensionedConstants"), dictionary(), false);
        controlDictPtr_ = new dictionary(dict);
    }

    return *controlDictPtr_;
}
*/


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        entry* ePtr = controlDict().lookupEntryPtr
        (
            subDictName, false, false
        );

        if (!ePtr || !ePtr->isDict())
        {
            cerr<< "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << controlDict().name().c_str()
                << std::endl << std::endl;

            ::exit(1);
        }
        subDictPtr = &ePtr->dict();
    }

    return *subDictPtr;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


Foam::dictionary& Foam::debug::tolerances()
{
    return switchSet("Tolerances", tolerancesPtr_);
}


int Foam::debug::debugSwitch(const char* name, const int defaultValue)
{
    return debugSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::infoSwitch(const char* name, const int defaultValue)
{
    return infoSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::optimisationSwitch(const char* name, const int defaultValue)
{
    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


float Foam::debug::floatOptimisationSwitch
(
    const char* name,
    const float defaultValue
)
{
    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


double Foam::debug::tolerances
(
    const char* name,
    const double defaultValue
)
{
    return tolerances().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


void Foam::debug::addDebugObject(const char* name, simpleRegIOobject* obj)
{
    simpleObjectRegistryEntry* ptr = debugObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        debugObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addInfoObject(const char* name, simpleRegIOobject* obj)
{
    simpleObjectRegistryEntry* ptr = infoObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        infoObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addOptimisationObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = optimisationObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        optimisationObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addDimensionSetObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = dimensionSetObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        dimensionSetObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addDimensionedConstantObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = dimensionedConstantObjects().lookupPtr
    (
        name
    );
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        dimensionedConstantObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


Foam::simpleObjectRegistry& Foam::debug::debugObjects()
{
    if (!debugObjectsPtr_)
    {
        debugObjectsPtr_ = new simpleObjectRegistry(1000);
    }

    return *debugObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::infoObjects()
{
    if (!infoObjectsPtr_)
    {
        infoObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *infoObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::optimisationObjects()
{
    if (!optimisationObjectsPtr_)
    {
        optimisationObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *optimisationObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionSetObjects()
{
    if (!dimensionSetObjectsPtr_)
    {
        dimensionSetObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *dimensionSetObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionedConstantObjects()
{
    if (!dimensionedConstantObjectsPtr_)
    {
        dimensionedConstantObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *dimensionedConstantObjectsPtr_;
}


// ************************************************************************* //
