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
    (c) 2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Table2.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::interpolation2DTable<Type>
Foam::Function2s::Table<Type>::readValues
(
    const word& tableName,
    const dictionary& dict
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>> values;
    Istream& is(dict.lookup(tableName));
    word entryType(is);
    is  >> values;
    return
        interpolation2DTable<Type>
        (
            values,
            interpolation2DTable<Type>::wordToBoundsHandling
            (
                dict.lookup("outOfBounds")
            ),
            word::null
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Table<Type>::Table
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction2<Type, Table<Type>>(name),
    values_(dict.found("file") ? (dict) : readValues(name, dict)),
    bounds_(dict.lookup("outOfBounds")),
    fileName_(dict.lookupOrDefault<word>("file", word::null))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function2s::Table<Type>::value
(
    scalar x,
    scalar y
) const
{
    return values_(x, y);
}


template<class Type>
void Foam::Function2s::Table<Type>::writeData(Ostream& os) const
{
    Function2<Type>::writeData(os);
    if (fileName_ != word::null)
    {
        os << token::END_STATEMENT << nl;
        os.writeEntry("file", fileName_);
    }
    else
    {
        os << values_;
        os << token::END_STATEMENT << nl;
    }
    os.writeEntry("outOfBounds", bounds_);
}


// ************************************************************************* //
