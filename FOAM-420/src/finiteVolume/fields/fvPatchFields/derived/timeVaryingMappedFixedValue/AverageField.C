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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/timeVaryingMappedFixedValue/AverageField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AverageField<Type>::AverageField(const label size)
:
    Field<Type>(size),
    average_(Zero)
{}


template<class Type>
Foam::AverageField<Type>::AverageField
(
    const Field<Type>& f,
    const Type& average
)
:
    Field<Type>(f),
    average_(average)
{}


template<class Type>
Foam::AverageField<Type>::AverageField(Istream& is)
:
    Field<Type>(is),
    average_(pTraits<Type>(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Type& Foam::AverageField<Type>::average() const
{
    return average_;
}


template<class Type>
Type&Foam::AverageField<Type>::average()
{
    return average_;
}


template<class Type>
bool Foam::AverageField<Type>::writeData(Ostream& os) const
{
    os  << static_cast<const Field<Type>&>(*this)
        << token::NL
        << average_;

    return os.good();
}


// ************************************************************************* //
