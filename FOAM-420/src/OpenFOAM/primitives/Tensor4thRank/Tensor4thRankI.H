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

#include "primitives/Tensor/Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<class Cmpt>
inline Tensor4thRank<Cmpt>::Tensor4thRank()
{}


// Construct given VectorSpace
template<class Cmpt>
inline Tensor4thRank<Cmpt>::Tensor4thRank
(
    const VectorSpace<Tensor4thRank<Cmpt>, Cmpt, 81>& vs
)
:
    VectorSpace<Tensor4thRank<Cmpt>, Cmpt, 81>(vs)
{}


// Construct given three Cmpts
template<class Cmpt>
inline Tensor4thRank<Cmpt>::Tensor4thRank
(
    const Cmpt xxxx, const Cmpt xxxy, const Cmpt xxxz, const Cmpt xxyx, const Cmpt xxyy, const Cmpt xxyz, const Cmpt xxzx, const Cmpt xxzy, const Cmpt xxzz,
    const Cmpt xyxx, const Cmpt xyxy, const Cmpt xyxz, const Cmpt xyyx, const Cmpt xyyy, const Cmpt xyyz, const Cmpt xyzx, const Cmpt xyzy, const Cmpt xyzz,
    const Cmpt xzxx, const Cmpt xzxy, const Cmpt xzxz, const Cmpt xzyx, const Cmpt xzyy, const Cmpt xzyz, const Cmpt xzzx, const Cmpt xzzy, const Cmpt xzzz,

    const Cmpt yxxx, const Cmpt yxxy, const Cmpt yxxz, const Cmpt yxyx, const Cmpt yxyy, const Cmpt yxyz, const Cmpt yxzx, const Cmpt yxzy, const Cmpt yxzz,
    const Cmpt yyxx, const Cmpt yyxy, const Cmpt yyxz, const Cmpt yyyx, const Cmpt yyyy, const Cmpt yyyz, const Cmpt yyzx, const Cmpt yyzy, const Cmpt yyzz,
    const Cmpt yzxx, const Cmpt yzxy, const Cmpt yzxz, const Cmpt yzyx, const Cmpt yzyy, const Cmpt yzyz, const Cmpt yzzx, const Cmpt yzzy, const Cmpt yzzz,

    const Cmpt zxxx, const Cmpt zxxy, const Cmpt zxxz, const Cmpt zxyx, const Cmpt zxyy, const Cmpt zxyz, const Cmpt zxzx, const Cmpt zxzy, const Cmpt zxzz,
    const Cmpt zyxx, const Cmpt zyxy, const Cmpt zyxz, const Cmpt zyyx, const Cmpt zyyy, const Cmpt zyyz, const Cmpt zyzx, const Cmpt zyzy, const Cmpt zyzz,
    const Cmpt zzxx, const Cmpt zzxy, const Cmpt zzxz, const Cmpt zzyx, const Cmpt zzyy, const Cmpt zzyz, const Cmpt zzzx, const Cmpt zzzy, const Cmpt zzzz
)
{
    this->v_[XXXX] = xxxx; this->v_[XXXY] = xxxy; this->v_[XXXZ] = xxxz; this->v_[XXYX] = xxyx; this->v_[XXYY] = xxyy; this->v_[XXYZ] = xxyz; this->v_[XXZX] = xxzx; this->v_[XXZY] = xxzy; this->v_[XXZZ] = xxzz;
    this->v_[XYXX] = xyxx; this->v_[XYXY] = xyxy; this->v_[XYXZ] = xyxz; this->v_[XYYX] = xyyx; this->v_[XYYY] = xyyy; this->v_[XYYZ] = xyyz; this->v_[XYZX] = xyzx; this->v_[XYZY] = xyzy; this->v_[XYZZ] = xyzz;
    this->v_[XZXX] = xzxx; this->v_[XZXY] = xzxy; this->v_[XZXZ] = xzxz; this->v_[XZYX] = xzyx; this->v_[XZYY] = xzyy; this->v_[XZYZ] = xzyz; this->v_[XZZX] = xzzx; this->v_[XZZY] = xzzy; this->v_[XZZZ] = xzzz;

    this->v_[YXXX] = yxxx; this->v_[YXXY] = yxxy; this->v_[YXXZ] = yxxz; this->v_[YXYX] = yxyx; this->v_[YXYY] = yxyy; this->v_[YXYZ] = yxyz; this->v_[YXZX] = yxzx; this->v_[YXZY] = yxzy; this->v_[YXZZ] = yxzz;
    this->v_[YYXX] = yyxx; this->v_[YYXY] = yyxy; this->v_[YYXZ] = yyxz; this->v_[YYYX] = yyyx; this->v_[YYYY] = yyyy; this->v_[YYYZ] = yyyz; this->v_[YYZX] = yyzx; this->v_[YYZY] = yyzy; this->v_[YYZZ] = yyzz;
    this->v_[YZXX] = yzxx; this->v_[YZXY] = yzxy; this->v_[YZXZ] = yzxz; this->v_[YZYX] = yzyx; this->v_[YZYY] = yzyy; this->v_[YZYZ] = yzyz; this->v_[YZZX] = yzzx; this->v_[YZZY] = yzzy; this->v_[YZZZ] = yzzz;

    this->v_[ZXXX] = zxxx; this->v_[ZXXY] = zxxy; this->v_[ZXXZ] = zxxz; this->v_[ZXYX] = zxyx; this->v_[ZXYY] = zxyy; this->v_[ZXYZ] = zxyz; this->v_[ZXZX] = zxzx; this->v_[ZXZY] = zxzy; this->v_[ZXZZ] = zxzz;
    this->v_[ZYXX] = zyxx; this->v_[ZYXY] = zyxy; this->v_[ZYXZ] = zyxz; this->v_[ZYYX] = zyyx; this->v_[ZYYY] = zyyy; this->v_[ZYYZ] = zyyz; this->v_[ZYZX] = zyzx; this->v_[ZYZY] = zyzy; this->v_[ZYZZ] = zyzz;
    this->v_[ZZXX] = zzxx; this->v_[ZZXY] = zzxy; this->v_[ZZXZ] = zzxz; this->v_[ZZYX] = zzyx; this->v_[ZZYY] = zzyy; this->v_[ZZYZ] = zzyz; this->v_[ZZZX] = zzzx; this->v_[ZZZY] = zzzy; this->v_[ZZZZ] = zzzz;
}


// Construct from Istream
template<class Cmpt>
inline Tensor4thRank<Cmpt>::Tensor4thRank(Istream& is)
:
    VectorSpace<Tensor4thRank<Cmpt>, Cmpt, 81>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
class outerProduct<Tensor4thRank<Cmpt>, Cmpt>
{
public:

    typedef Tensor4thRank<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, Tensor4thRank<Cmpt>>
{
public:

    typedef Tensor4thRank<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Tensor4thRank<Cmpt>, Vector<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, Tensor4thRank<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class scalarProduct<Tensor4thRank<Cmpt>, Tensor<Cmpt>>
{
public:

  typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class scalarProduct<Tensor<Cmpt>, Tensor4thRank<Cmpt>>
{
public:

  typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
