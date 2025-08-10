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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2019-2021 OpenCFD Ltd.

Class
    Foam::simpleObjectRegistry

Description
    Object registry for simpleRegIOobject. Maintains ordering.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef simpleObjectRegistry_H
#define simpleObjectRegistry_H

#include "containers/Dictionaries/Dictionary/Dictionary.H"
#include "global/debug/simpleRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                   Class simpleObjectRegistryEntry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistryEntry
:
    public Dictionary<simpleObjectRegistryEntry>::link,
    public List<simpleRegIOobject*>
{
public:

    simpleObjectRegistryEntry(const List<simpleRegIOobject*>& data)
    :
        List<simpleRegIOobject*>(data)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class simpleObjectRegistry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistry
:
    public Dictionary<simpleObjectRegistryEntry>
{
public:

    // Constructors

        //- Construct given initial table size
        simpleObjectRegistry(const label nIoObjects = 128)
        :
            Dictionary<simpleObjectRegistryEntry>(nIoObjects)
        {}


    // Member Functions

        //- Set named value, but also handle embedded 'name=value' syntax
        //  Treats 'name=' and 'name' identically.
        //
        //  Will also accept a floating point value, but the called objects
        //  need to handle that.
        void setNamedValue(std::string name, int val, bool report=false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
