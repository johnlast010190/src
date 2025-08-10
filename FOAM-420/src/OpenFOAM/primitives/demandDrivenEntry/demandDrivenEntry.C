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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/demandDrivenEntry/demandDrivenEntry.H"

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const Type& value
)
:
    dict_(dict),
    keyword_("unknown-keyword"),
    value_(value),
    stored_(true)
{}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const word& keyword
)
:
    dict_(dict),
    keyword_(keyword),
    value_(Zero),
    stored_(false)
{}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry
(
    const dictionary& dict,
    const word& keyword,
    const Type& defaultValue,
    const bool readIfPresent
)
:
    dict_(dict),
    keyword_(keyword),
    value_(defaultValue),
    stored_(true)
{
    if (readIfPresent)
    {
        dict_.readIfPresent<Type>(keyword, value_);
    }
}


template<class Type>
Foam::demandDrivenEntry<Type>::demandDrivenEntry(const demandDrivenEntry& dde)
:
    dict_(dde.dict_),
    keyword_(dde.keyword_),
    value_(dde.value_),
    stored_(dde.stored_)
{}


// ************************************************************************* //
