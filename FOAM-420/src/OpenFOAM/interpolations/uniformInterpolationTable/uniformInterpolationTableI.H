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
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

template<class Type>
Foam::scalar Foam::uniformInterpolationTable<Type>::x0() const
{
    return x0_;
}


template<class Type>
Foam::scalar Foam::uniformInterpolationTable<Type>::dx() const
{
    return dx_;
}


template<class Type>
const Foam::Switch& Foam::uniformInterpolationTable<Type>::log10() const
{
    return log10_;
}


template<class Type>
const Foam::Switch& Foam::uniformInterpolationTable<Type>::bound() const
{
    return bound_;
}


template<class Type>
Foam::scalar& Foam::uniformInterpolationTable<Type>::x0()
{
    return x0_;
}


template<class Type>
Foam::scalar& Foam::uniformInterpolationTable<Type>::dx()
{
    return dx_;
}


template<class Type>
Foam::Switch& Foam::uniformInterpolationTable<Type>::log10()
{
    return log10_;
}


template<class Type>
Foam::Switch& Foam::uniformInterpolationTable<Type>::bound()
{
    return bound_;
}


template<class Type>
Foam::scalar Foam::uniformInterpolationTable<Type>::xMin() const
{
    return x0_;
}


template<class Type>
Foam::scalar Foam::uniformInterpolationTable<Type>::xMax() const
{
    return x0_ + dx_*(size() - 1);
}


// ************************************************************************* //
