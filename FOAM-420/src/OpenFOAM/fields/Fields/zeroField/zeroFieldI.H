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

\*---------------------------------------------------------------------------*/

#include "fields/Fields/zeroField/zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::zeroField::operator[](const label) const
{
    return scalar(0);
}


inline Foam::zeroField Foam::zeroField::field() const
{
    return zeroField();
}


inline Foam::zeroField Foam::zeroField::operator()() const
{
    return zeroField();
}


inline Foam::zeroField Foam::zeroField::operator-() const
{
    return zeroField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const zeroField operator+(const zeroField&, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline const Type& operator+(const Type& t, const zeroField&)
{
    return t;
}

template<class Type>
inline const Type& operator+(const zeroField&, const Type& t)
{
    return t;
}

inline const zeroField operator-(const zeroField&, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline const Type& operator-(const Type& t, const zeroField&)
{
    return t;
}

template<class Type>
inline Type operator-(const zeroField&, const Type& t)
{
    return -t;
}

template<class Type>
inline zeroField operator*(const Type& t, const zeroField&)
{
    return zeroField();
}

template<class Type>
inline zeroField operator*(const zeroField&, const Type& t)
{
    return zeroField();
}

template<class Type>
inline zeroField operator/(const zeroField&, const Type& t)
{
    return zeroField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
