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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "removeRegisteredObject/removeRegisteredObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(removeRegisteredObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        removeRegisteredObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::removeRegisteredObject::removeRegisteredObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::removeRegisteredObject::~removeRegisteredObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::removeRegisteredObject::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    dict.lookup("objects") >> objectNames_;

    return true;
}


bool Foam::functionObjects::removeRegisteredObject::execute()
{
    forAll(objectNames_, i)
    {
        if (foundObject<regIOobject>(objectNames_[i]))
        {
            const regIOobject& obj =
                lookupObject<regIOobject>(objectNames_[i]);

            if (obj.ownedByRegistry())
            {
                Log << type() << " " << name() << " output:" << nl
                    << "    removing object " << obj.name() << nl
                    << endl;

                const_cast<regIOobject&>(obj).release();
                delete &obj;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::removeRegisteredObject::write()
{
    return true;
}


// ************************************************************************* //
