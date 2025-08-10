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
    (c) 2019 Esi Ltd.

InNamespace
    Foam

Description
    Traits classes for inner and outer products of primitives.

\*---------------------------------------------------------------------------*/

#ifndef products_H
#define products_H

#include "primitives/pTraits/pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Abstract template class to provide the form resulting from
//  the inner-product of two forms
template<class Cmpt, class Form1, class Form2>
class typeOfInnerProduct
{};

//- Abstract template class to provide the form resulting from
//  the outer-product of two forms
template<class Cmpt, class Form1, class Form2>
class typeOfOuterProduct
{};

//- Abstract template class to provide the transpose form of a form
template<class Cmpt, class Form>
class typeOfTranspose
{};


template<class Cmpt, int rank>
class typeOfRank
{};


template<class Cmpt, int rank>
class symmTypeOfRank
{};


template<class arg1, class arg2>
class typeOfSum
{
public:

    typedef arg1 type;
};


template<class arg1, class arg2>
class outerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank)
    >::type type;
};


template<class arg1, class arg2>
class crossProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg2>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank) - 1
    >::type type;
};

template<class arg1, class arg2>
class innerProduct
{
public:

    typedef typename typeOfRank
    <
        typename pTraits<arg1>::cmptType,
        int(pTraits<arg1>::rank) + int(pTraits<arg2>::rank) - 2
    >::type type;
};

template<class arg1, class arg2>
class scalarProduct
{
public:

    typedef typename pTraits<arg1>::cmptType type;
};


template<class arg1, int arg2>
class powProduct
{
public:

    typedef typename symmTypeOfRank
    <
        typename pTraits<arg1>::cmptType,
        arg2*int(pTraits<arg1>::rank)
    >::type type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
