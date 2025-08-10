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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    IOLongList

Description
    An IOLongList of a given type is a List of that type which supports
    automated input and output.

SourceFiles
    IOLongList.C

\*---------------------------------------------------------------------------*/

#ifndef IOLongList_H
#define IOLongList_H

#include "utilities/containers/LongList/LongList.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOLongList Declaration
\*---------------------------------------------------------------------------*/

template<class T, label Offset = 19>
class IOLongList
:
    public regIOobject,
    public LongList<T, Offset>
{

public:

    //- Runtime type information
    TypeName("LongList");


    // Constructors

        //- Construct from IOobject
        IOLongList(const IOobject&);

        //- Construct from IOobject and size of IOLongList
        IOLongList(const IOobject&, const label);

        //- Construct from IOobject and a List
        IOLongList(const IOobject&, const LongList<T, Offset>&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOLongList<T, Offset>&);

        void operator=(const LongList<T, Offset>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "utilities/containers/IOLongList/IOLongList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
