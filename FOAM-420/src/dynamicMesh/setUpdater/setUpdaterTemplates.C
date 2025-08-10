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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "setUpdater/setUpdater.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::setUpdater::updateSets(const mapPolyMesh& morphMap) const
{
    //
    // Update all sets in memory.
    //

    HashTable<const Type*> memSets =
        morphMap.mesh().objectRegistry::lookupClass<Type>();

    forAllIter(typename HashTable<const Type*>, memSets, iter)
    {
        Type& set = const_cast<Type&>(*iter());

        if (debug)
        {
            Pout<< "Set:" << set.name() << " size:" << set.size()
                << " updated in memory" << endl;
        }

        set.updateMesh(morphMap);

        // Write or not? Debatable.
        set.write();
    }


    //
    // Update all sets on disk
    //

    // Get last valid mesh (discard points-only change)
    IOobjectList Objects
    (
        morphMap.mesh().time(),
        morphMap.mesh().facesInstance(),
        "polyMesh/sets"
    );

    IOobjectList fileSets(Objects.lookupClass(Type::typeName));

    forAllConstIter(IOobjectList, fileSets, iter)
    {
        if (!memSets.found(iter.key()))
        {
            // Not in memory. Load it.
            Type set(*iter());

            if (debug)
            {
                Pout<< "Set:" << set.name() << " size:" << set.size()
                    << " updated on disk" << endl;
            }

            set.updateMesh(morphMap);

            set.write();
        }
        else
        {
            if (debug)
            {
                Pout<< "Set:" << iter.key() << " already updated from memory"
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
