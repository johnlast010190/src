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

Class
    Foam::boundaryInfo

Description
    Class to interrogate the polyMesh/boundary file to provide mesh patching
    information, without the need to read the mesh.

\*---------------------------------------------------------------------------*/

#ifndef boundaryInfo_H
#define boundaryInfo_H

#include "primitives/bools/lists/boolList.H"
#include "primitives/strings/lists/wordList.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "db/IOobjects/IOPtrList/IOPtrList.H"
#include "db/dictionary/entry/entry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                        Class boundaryInfo Declaration
\*---------------------------------------------------------------------------*/

class boundaryInfo
{
    // Private data

        //- Boundary dictionary
        IOPtrList<entry> boundaryDict_;

        //- Patch names
        wordList names_;

        //- Patch types
        wordList types_;

        //- Constraint flag
        boolList constraint_;

        //- Groups per patch
        List<wordList> groups_;

        //- Set of all group names
        wordHashSet allGroupNames_;


    // Private member functions

        //- Read the boundary dict
        IOPtrList<entry> readBoundaryDict
        (
            const Time& runTime,
            const word& regionName
        ) const;


public:

    //- Constructor
    boundaryInfo(const Time& runTime, const word& regionName);


    // Public member functions

        //- Patch names
        const wordList& names() const;

        //- Patch types
        const wordList& types() const;

        //- Constraint flag
        const boolList& constraint() const;

        //- Groups
        const List<wordList>& groups() const;

        //- Set of all group names
        const wordHashSet& allGroupNames() const;

        //- Set the patch type based on the condition
        void setType(const label patchI, const word& condition);

        //- Write the boundary dictionary
        void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
