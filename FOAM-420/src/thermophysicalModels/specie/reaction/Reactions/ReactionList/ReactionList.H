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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::ReactionList

Description
    List of templated reactions

SourceFiles
    ReactionList.C

\*---------------------------------------------------------------------------*/

#ifndef ReactionList_H
#define ReactionList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "containers/LinkedLists/user/SLPtrList.H"
#include "speciesTable/speciesTable.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "reaction/Reactions/Reaction/Reaction.H"
#include "primitives/strings/fileName/fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ReactionList Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class ReactionList
:
    public SLPtrList<Reaction<ThermoType>>
{
    // Private data

        //- Reference to the table of species
        const speciesTable& species_;

        //- Reference to the thermo database
        const HashPtrTable<ThermoType>& thermoDb_;

        //- The dictionary used for construction
        const dictionary dict_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ReactionList&);


public:

    // Constructors

        //- Construct null
        ReactionList
        (
            const speciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase
        );

        //- Construct from dictionary
        ReactionList
        (
            const speciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct copy
        ReactionList(const ReactionList<ThermoType>& reactions);


    //- Destructor
    ~ReactionList();


    // Public Member Functions

        //- Read reactions from dictionary
        bool readReactionDict();

        //- Write
        void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reaction/Reactions/ReactionList/ReactionList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
