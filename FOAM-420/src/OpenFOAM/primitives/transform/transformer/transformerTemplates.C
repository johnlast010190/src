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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::transformer::transform(const Type& x) const
{
    if (transforms())
    {
        return Foam::transform(T(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
void Foam::transformer::transform
(
    Field<Type>& res,
    const Field<Type>& fld
) const
{
    if (transforms())
    {
        return Foam::transform(res, T(), fld);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::transform
(
    const Field<Type>& fld
) const
{
    if (transforms())
    {
        return Foam::transform(T(), fld);
    }
    else
    {
        return fld;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::transform
(
    const tmp<Field<Type>>& tfld
) const
{
    if (transforms())
    {
        return Foam::transform(T(), tfld);
    }
    else
    {
        return tfld;
    }
}


template<class Type, template<class> class Container>
void Foam::transformer::transformList(Container<Type>& l) const
{
    if (transforms())
    {
        forAllIter(typename Container<Type>, l, iter)
        {
            *iter = Foam::transform(T(), *iter);
        }
    }
}


template<class Type>
Type Foam::transformer::invTransform(const Type& x) const
{
    if (transforms())
    {
        return Foam::transform(invT(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
void Foam::transformer::invTransform
(
    Field<Type>& res,
    const Field<Type>& fld
) const
{
    if (transforms())
    {
        Foam::transform(res, invT(), fld);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::invTransform
(
    const Field<Type>& fld
) const
{
    if (transforms())
    {
        return Foam::transform(invT(), fld);
    }
    else
    {
        return fld;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::invTransform
(
    const tmp<Field<Type>>& tfld
) const
{
    if (transforms())
    {
        return Foam::transform(invT(), tfld);
    }
    else
    {
        return tfld;
    }
}


template<class Type, template<class> class Container>
void Foam::transformer::invTransformList(Container<Type>& l) const
{
    if (transforms())
    {
        tensor invT = this->invT();

        forAllIter(typename Container<Type>, l, iter)
        {
            *iter = Foam::transform(invT, *iter);
        }
    }
}


// ************************************************************************* //
