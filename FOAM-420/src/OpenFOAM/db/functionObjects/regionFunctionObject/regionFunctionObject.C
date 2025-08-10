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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"
#include "db/solutionRegistry/solutionRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(regionFunctionObject, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::whichRegistry
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    HashTable<const solutionInstanceRegistry*> instanceRegistries =
        obr.lookupClass<solutionInstanceRegistry>();
    word regionName = dict.lookupOrDefault("region", polyMesh::defaultRegion);
    if (instanceRegistries.size())
    {
        word meshName;
        if (dict.found("instance"))
        {
            word instanceName = dict.lookup<word>("instance");
            const auto& solReg =
                obr.lookupObject<solutionInstanceRegistry>(instanceName);
            meshName = solReg.whichMesh(regionName);
        }
        else
        {
            // Use the active instance
            forAllIters(instanceRegistries, iter)
            {
                if (iter()->isActive())
                {
                    meshName = iter()->whichMesh(regionName);
                    break;
                }
            }
        }
        if (obr.foundObject<objectRegistry>(meshName))
        {
            const objectRegistry& meshDb =
                obr.lookupObject<objectRegistry>(meshName);
            return meshDb;
        }
        else
        {
            end();
            return obr;
        }
    }
    else
    {
        return obr.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        );
    }
}

const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::whichSubRegistry
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    HashTable<const solutionInstanceRegistry*> instanceRegistries =
        obr.lookupClass<solutionInstanceRegistry>();
    word subName = dict.lookupOrDefault("region", polyMesh::defaultRegion);
    if (instanceRegistries.size() && subName != polyMesh::defaultRegion)
    {
        word meshName = meshObr_.name();
        return
            meshObr_.lookupObject<solutionRegistry>
            (
                (subName == meshName) ? "" : subName
            ).registry();
    }
    else
    {
        return meshObr_;
    }
}


const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::obr() const
{
    return obr_;
}


bool Foam::functionObjects::regionFunctionObject::writeObject
(
    const word& fieldName
)
{
    const regIOobject* objPtr =
        this->lookupObjectPtr<regIOobject>(fieldName);

    if (objPtr)
    {
        Log << "    functionObjects::" << type() << " " << name()
            << " writing field: " << objPtr->name() << endl;

        objPtr->write();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::regionFunctionObject::clearObject
(
    const word& fieldName
)
{
    regIOobject* objPtr = lookupObjectRefPtr<regIOobject>(fieldName);
    if (objPtr)
    {
        if (objPtr->ownedByRegistry())
        {
            return objPtr->checkOut();
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    meshObr_(whichRegistry(runTime, dict)),
    obr_(whichSubRegistry(runTime, dict))
{}


Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    stateFunctionObject(name, obr.time()),
    meshObr_(obr),
    obr_(whichSubRegistry(obr.time(), dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::~regionFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionFunctionObject::read(const dictionary& dict)
{
    return stateFunctionObject::read(dict);
}


// ************************************************************************* //
