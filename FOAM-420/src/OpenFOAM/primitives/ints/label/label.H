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
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

Typedef
    Foam::label

Description
    A label is an int32_t or int64_t as specified by the pre-processor macro
    FOAM_LABEL_SIZE.

    A readLabel function is defined so that label can be constructed from
    Istream.

\*---------------------------------------------------------------------------*/

#ifndef label_H
#define label_H

#include "primitives/ints/int/int.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INT_ADD_SIZE(x,s,y) x ## s ## y
#define INT_ADD_DEF_SIZE(x,s,y) INT_ADD_SIZE(x,s,y)
#define INT_SIZE(x,y) INT_ADD_DEF_SIZE(x,FOAM_LABEL_SIZE,y)

#if FOAM_LABEL_SIZE != 32 && FOAM_LABEL_SIZE != 64
    #error "label.H: FOAM_LABEL_SIZE must be set to either 32 or 64"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef INT_SIZE(int, _t) label;

static const label labelMin = INT_SIZE(INT, _MIN);
static const label labelMax = INT_SIZE(INT, _MAX);

inline label readLabel(Istream& is)
{
    return INT_SIZE(readInt,) (is);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one label to the power of another
label pow(label a, label b);

//- Evaluate n! : 0 < n <= 12
label factorial(label n);


inline label& setComponent(label& l, const direction)
{
    return l;
}

inline label component(const label l, const direction)
{
    return l;
}

inline label sign101(const label s) //- dcombest 8/15/2012
{
    return (s >= 0)-(s <= 0);
}

/*inline label sign101(const label s, const label w)//- dcombest 10/8/2012
{
    return s/pow((s*s + w*w + 1e-20),0.5);
}*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/ints/label/labelSpecific.H"

#undef INT_ADD_SIZE
#undef INT_ADD_DEF_SIZE
#undef INT_SIZE

#endif

// ************************************************************************* //
