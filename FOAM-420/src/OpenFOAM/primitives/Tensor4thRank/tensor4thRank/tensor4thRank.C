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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/Tensor4thRank/tensor4thRank/tensor4thRank.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const tensor4thRank::typeName = "tensor4thRank";

template<>
const char* tensor4thRank::componentNames[] =
{
    "xxxx", "xxxy", "xxxz", "xxyx", "xxyy", "xxyz", "xxzx", "xxzy", "xxzz",
    "xyxx", "xyxy", "xyxz", "xyyx", "xyyy", "xyyz", "xyzx", "xyzy", "xyzz",
    "xzxx", "xzxy", "xzxz", "xzyx", "xzyy", "xzyz", "xzzx", "xzzy", "xzzz",

    "yxxx", "yxxy", "yxxz", "yxyx", "yxyy", "yxyz", "yxzx", "yxzy", "yxzz",
    "yyxx", "yyxy", "yyxz", "yyyx", "yyyy", "yyyz", "yyzx", "yyzy", "yyzz",
    "yzxx", "yzxy", "yzxz", "yzyx", "yzyy", "yzyz", "yzzx", "yzzy", "yzzz",

    "zxxx", "zxxy", "zxxz", "zxyx", "zxyy", "zxyz", "zxzx", "zxzy", "zxzz"
    "zyxx", "zyxy", "zyxz", "zyyx", "zyyy", "zyyz", "zyzx", "zyzy", "zyzz"
    "zzxx", "zzxy", "zzxz", "zzyx", "zzyy", "zzyz", "zzzx", "zzzy", "zzzz"
};

template<>
const tensor4thRank tensor4thRank::zero
(
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0
);

template<>
const tensor4thRank tensor4thRank::one
(
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,

    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,

    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1
);

template<>
const tensor4thRank tensor4thRank::max
(
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,

    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,

    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT, VGREAT
);

template<>
const tensor4thRank tensor4thRank::min
(
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,

    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,

    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
