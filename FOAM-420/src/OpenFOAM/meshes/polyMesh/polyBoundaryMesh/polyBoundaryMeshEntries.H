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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::polyBoundaryMeshEntries

Description
    Foam::polyBoundaryMeshEntries

SourceFiles
    polyBoundaryMeshEntries.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMeshEntries_H
#define polyBoundaryMeshEntries_H

#include "db/regIOobject/regIOobject.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/dictionary/entry/entry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyBoundaryMeshEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        explicit polyBoundaryMeshEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(typeName))
        {
            close();
        }


   // Member functions

        bool writeData(Ostream&) const
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
