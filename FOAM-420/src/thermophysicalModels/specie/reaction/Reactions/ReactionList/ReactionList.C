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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "reaction/Reactions/ReactionList/ReactionList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "containers/LinkedLists/user/SLPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb
)
:
    SLPtrList<Reaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const dictionary& dict
)
:
    SLPtrList<Reaction<ThermoType>>(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dict)
{
    readReactionDict();
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList(const ReactionList& reactions)
:
    SLPtrList<Reaction<ThermoType>>(reactions),
    species_(reactions.species_),
    thermoDb_(reactions.thermoDb_),
    dict_(reactions.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::~ReactionList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::ReactionList<ThermoType>::readReactionDict()
{
    if (dict_.found("reactions"))
    {
        const dictionary& reactions(dict_.subDict("reactions"));

        forAllConstIter(dictionary, reactions, iter)
        {
            const word reactionName = iter().keyword();

            this->append
            (
                Reaction<ThermoType>::New
                (
                    species_,
                    thermoDb_,
                    reactions.subDict(reactionName)
                ).ptr()
            );
        }
    }

    return true;
}


template<class ThermoType>
void Foam::ReactionList<ThermoType>::write(Ostream& os) const
{
    os.beginBlock("reactions");
    forAllConstIter(typename SLPtrList<Reaction<ThermoType>>, *this, iter)
    {
        const Reaction<ThermoType>& r = iter();
        os.beginBlock(r.name());
        os.writeEntry("type", r.type());
        r.write(os);
        os.endBlock();
    }
    os.endBlock();
}


// ************************************************************************* //
