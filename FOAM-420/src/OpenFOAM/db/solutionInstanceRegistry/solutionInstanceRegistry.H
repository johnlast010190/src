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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef solutionInstanceRegistry_H
#define solutionInstanceRegistry_H

#include "db/Time/Time.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/lists/hashedWordList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionInstanceRegistry
:
    public objectRegistry
{

    //- Disallow default bitwise copy construct
    solutionInstanceRegistry(const solutionInstanceRegistry&);

    //- Disallow default bitwise assignment
    void operator=(const solutionInstanceRegistry&);


protected:

    dictionary solutionMeshDict_;

    DynamicList<word> meshNames_;

    List<wordList> solutionRegions_;

    wordList regionNames_;

    mutable bool active_;

    //- Flag to tell which instances are selected for the scheduler
    mutable bool selected_;

    label getTotalNumberOfRegions();

    void generateMeshList();

    label whichRegionIndex(const word& regionName) const;


public:

    // Constructors
        solutionInstanceRegistry
        (
            const objectRegistry& obr,
            const List<List<word>>& instanceRegions,
            const dictionary& dict,
            const word& name
        );

    ~solutionInstanceRegistry();

    // Member Functions
        const word whichMesh(const word& regionName) const;

        inline const bool& isActive() const
        {
            return active_;
        }

        inline const bool& isSelected() const
        {
            return selected_;
        }

        inline const wordList& regionNames() const
        {
            return regionNames_;
        }

        inline const wordList& meshNames() const
        {
            return meshNames_;
        }

        void unselectInstance()
        {
            selected_ = false;
        }
};

}

#endif
