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

InClass
    Foam::makeSolidReaction

Description
    Macros for instantiating reactions for solid

\*---------------------------------------------------------------------------*/

#ifndef makeSolidReactionThermo_H
#define makeSolidReactionThermo_H

#include "reaction/Reactions/solidReaction/solidReaction.H"
#include "reaction/Reactions/IrreversibleReaction/IrreversibleReaction.H"
#include "reaction/Reactions/Reaction/Reaction.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidReaction(ReactionType, Thermo, ReactionRate)                  \
                                                                               \
    typedef solidReaction<Thermo> solidReaction##Thermo;                       \
                                                                               \
    typedef Reaction<Thermo> Reaction##Thermo;                                 \
                                                                               \
    typedef ReactionType<solidReaction, Thermo, ReactionRate>                  \
        ReactionType##Thermo##ReactionRate;                                    \
                                                                               \
    defineTemplateRunTimeSelectionTable(Reaction##Thermo, dictionary);         \
                                                                               \
    defineTemplateTypeNameAndDebug(solidReaction##Thermo, 0);                  \
    defineTemplateTypeNameAndDebug(Reaction##Thermo, 0);                       \
                                                                               \
    template<>                                                                 \
    const word ReactionType##Thermo##ReactionRate::typeName                    \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type()                                                   \
      + solidReaction##Thermo::typeName_()                                     \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##Thermo,                                                      \
        ReactionType##Thermo##ReactionRate,                                    \
        dictionary                                                             \
    );

#define makeSolidIRReactions(Thermo, ReactionRate)                             \
                                                                               \
    makeSolidReaction(IrreversibleReaction, Thermo, ReactionRate)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
