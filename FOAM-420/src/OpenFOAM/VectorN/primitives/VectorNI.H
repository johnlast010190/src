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
    (c) 2016 Esi Ltd.
    (c) 2010 Ivor Clifford

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Cmpt, int length>
class TensorN;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Cmpt, int length>
const VectorN<Cmpt, length> VectorN<Cmpt, length>::zero(0);

template <class Cmpt, int length>
const VectorN<Cmpt, length> VectorN<Cmpt, length>::one(1);

template <class Cmpt, int length>
const VectorN<Cmpt, length> VectorN<Cmpt, length>::max(VGREAT);

template <class Cmpt, int length>
const VectorN<Cmpt, length> VectorN<Cmpt, length>::min(-VGREAT);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template <class Cmpt, int length>
inline VectorN<Cmpt, length>::VectorN()
{}

template <class Cmpt, int length>
inline VectorN<Cmpt, length>::VectorN(const Foam::zero)
:
   // VectorN::vsType(Zero)
    VectorSpace<VectorN<Cmpt, length>, Cmpt, length>(Zero)
{}

// Construct given VectorSpace
template <class Cmpt, int length>
inline VectorN<Cmpt, length>::VectorN
(
    const VectorSpace<VectorN<Cmpt, length>, Cmpt, length>& vs
)
:
    VectorSpace<VectorN<Cmpt, length>, Cmpt, length>(vs)
    //VectorN::vsType(vs)
{}


// Construct given Cmpts
template <class Cmpt, int length>
inline VectorN<Cmpt, length>::VectorN(const Cmpt& vx)
{
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::eqOpS(*this, vx, eqOp<Cmpt>());
}


// Construct from Istream
template <class Cmpt, int length>
inline VectorN<Cmpt, length>::VectorN(Istream& is)
:
    VectorSpace<VectorN<Cmpt, length>, Cmpt, length>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Cmpt, int length>
inline const Cmpt& VectorN<Cmpt, length>::operator()(const direction i) const
{
    return this->operator[](i);
}


template <class Cmpt, int length>
inline Cmpt& VectorN<Cmpt, length>::operator()(const direction i)
{
    return this->operator[](i);
}

//- Multiply components of VectorN by VectorN
template<class Cmpt, int length>
inline VectorN<Cmpt, length> VectorN<Cmpt, length>::cmptMultiply
(
    const VectorN<Cmpt, length>& v
)
{
    VectorN<Cmpt, length> res;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::op
        (res, *this, v, multiplyOp<Cmpt>());
    return res;
}


//- Transform the spherical tensor
//- The components are assumed to be individual scalars
//- i.e. transform has no effect
template<class Cmpt, int length>
inline VectorN<Cmpt, length> transform
(
    const tensor& tt,
    const VectorN<Cmpt, length>& v
)
{
    return v;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt, int length>
class outerProduct<Cmpt, VectorN<Cmpt, length>>
{
public:

    typedef VectorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class outerProduct<VectorN<Cmpt, length>, Cmpt>
{
public:

    typedef VectorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class outerProduct<VectorN<Cmpt, length>, VectorN<Cmpt, length>>
{
public:

    typedef TensorN<Cmpt, length> type;
};


//- Inner-product between a vector and a vector
template <class Cmpt, int length>
inline Cmpt operator&
(
    const VectorN<Cmpt, length>& v1,
    const VectorN<Cmpt, length>& v2
)
{
    Cmpt res = pTraits<Cmpt>::zero;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::SopEqOpVV
        (res, v1, v2, plusEqOp<Cmpt>(), multiplyOp<Cmpt>());
    return res;
}


//- Componentwise division of scalar by VectorN
template<class Cmpt, int length>
inline VectorN<Cmpt, length>
operator/(const scalar s, const VectorN<Cmpt, length>& v)
{
    VectorN<Cmpt, length> res;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::opSV
        (res, s, v, divideOp3<Cmpt, scalar, Cmpt>());
    return res;
}


//- Multiply components of VectorN by VectorN
template<class Cmpt, int length>
inline VectorN<Cmpt, length> cmptMultiply
(
    const VectorN<Cmpt, length>& v1,
    const VectorN<Cmpt, length>& v2
)
{
    VectorN<Cmpt, length> res;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::op
        (res, v1, v2, multiplyOp<Cmpt>());
    return res;
}

//- Return the component sum
// template <class Cmpt, int length>
// inline Cmpt sum(const VectorN<Cmpt, length>& v)
// {
//     Cmpt result = pTraits<Cmpt>::zero;
//     for (register label i=0; i<VectorN<Cmpt, length>::nComponents; i++)
//     {
//         result += v[i];
//     }
//     return result;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
