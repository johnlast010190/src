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
    (c) 2004-6 H. Jasak All rights reserved
    (c) 2019 Esi Ltd.

Class
    blockCoeff

Description
    Block coefficient combines a scalar, linear and square coefficient
    for different levels of coupling

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef BlockCoeff_H
#define BlockCoeff_H

#include "primitives/BlockCoeff/blockCoeffBase.H"
#include "fields/expandContract/expandTensor.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockCoeff;

template<class Type>
Ostream& operator<<(Ostream&, const BlockCoeff<Type>&);

/*---------------------------------------------------------------------------*\
                          Class BlockCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockCoeff
:
    public blockCoeffBase
{
public:

    // Public data types

        //- Component type
        typedef Type xType;
        typedef Field<xType> xTypeField;

        //- Coefficient type
        typedef typename pTraits<Type>::cmptType scalarType;
        typedef Type linearType;
        typedef typename outerProduct<Type, Type>::type squareType;

        //- Field type
        typedef Field<scalarType> scalarTypeField;
        typedef Field<linearType> linearTypeField;
        typedef Field<squareType> squareTypeField;


    //- Multiplication trait
    class multiply
    {
    public:
        // Coefficient times type multiplication

            Type operator()(const scalarType& c, const Type& x) const
            {
                return c*x;
            }

            Type operator()(const linearType& c, const Type& x) const
            {
                return cmptMultiply(c, x);
            }

            Type operator()(const squareType& c, const Type& x) const
            {
                return (c & x);
            }

            Type operator()(const BlockCoeff<Type>& c, const Type& x) const
            {
                if (c.scalarCoeffPtr_)
                {
                    return operator()(*c.scalarCoeffPtr_, x);
                }
                else if (c.linearCoeffPtr_)
                {
                    return operator()(*c.linearCoeffPtr_, x);
                }
                else if (c.squareCoeffPtr_)
                {
                    return operator()(*c.squareCoeffPtr_, x);
                }
                else
                {
                    return pTraits<Type>::zero;
                }
            }

            scalar operator()(const direction row, const scalarType& c, const Type& x) const
            {
                return c*x;
            }

            scalar operator()(const direction row, const linearType& c, const Type& x) const
            {
                return c(row)*x(row);
            }

            scalar operator()(const direction row, const squareType& c, const Type& x) const
            {
                scalar result(0);
                const direction nCmpts = pTraits<Type>::nComponents;
                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    result += c(row, cmptI)*x(cmptI);
                }
                return result;
            }

            scalar operator()(const direction row, const BlockCoeff<Type>& c, const Type& x) const
            {
                if (c.scalarCoeffPtr_)
                {
                    return operator()(row, *c.scalarCoeffPtr_, x);
                }
                else if (c.linearCoeffPtr_)
                {
                    return operator()(row, *c.linearCoeffPtr_, x);
                }
                else if (c.squareCoeffPtr_)
                {
                    return operator()(row, *c.squareCoeffPtr_, x);
                }
                else
                {
                    return pTraits<Type>::zero;
                }
            }

            scalarType activeTypeMultiply
            (
                const scalarType& a,
                const scalarType& b
            ) const
            {
                return a*b;
            }

            linearType activeTypeMultiply
            (
                const linearType& a,
                const linearType& b
            ) const
            {
                return cmptMultiply(a, b);
            }

            squareType activeTypeMultiply
            (
                const squareType& a,
                const squareType& b
            ) const
            {
                return (a & b);
            }


        // Transpose functions

            scalarType transpose(const scalarType& c) const
            {
                return c;
            }

            linearType transpose(const linearType& c) const
            {
                return c;
            }

            squareType transpose(const squareType& c) const
            {
                return c.T();
            }


        // Inverse functions

            scalarType inverse(const scalarType& c) const
            {
                return 1.0/c;
            }

            linearType inverse(const linearType& c) const
            {
                return cmptDivide(pTraits<linearType>::one, c);
            }

            squareType inverse(const squareType& c) const
            {
                return inv(c);
            }

            scalarType diagMultiply(const scalarType& a, const scalar& c) const
            {
                return a*c;
            }

            linearType diagMultiply(const linearType& a, const scalar& c) const
            {
                return a*c;
            }

            squareType diagMultiply(const squareType& a,  const scalar& c) const
            {
                squareType result = a;
                const direction nCmpts = pTraits<Type>::nComponents;
                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    result(cmptI, cmptI) = a(cmptI, cmptI)*c;
                }
                return result;
            }


        // Triple product of coefficients

            scalarType tripleProduct
            (
                const scalarType& a,
                const scalarType& b,
                const scalarType& c
            ) const
            {
                return a*c/b;
            }

            linearType tripleProduct
            (
                const scalarType& a,
                const linearType& b,
                const scalarType& c
            ) const
            {
                return a*c*inverse(b);
            }

            linearType tripleProduct
            (
                const linearType& a,
                const linearType& b,
                const linearType& c
            ) const
            {
                return cmptDivide(cmptMultiply(a, c), b);
            }

            squareType tripleProduct
            (
                const scalarType& a,
                const squareType& b,
                const scalarType& c
            ) const
            {
                return a*c*inv(b);
            }

            squareType tripleProduct
            (
                const linearType& a,
                const squareType& b,
                const linearType& c
            ) const
            {
                squareType result;
                linearType sac = cmptMultiply(a, c);

                expandLinear(result, sac);
                return result & inv(b);
            }

            squareType tripleProduct
            (
                const squareType& a,
                const squareType& b,
                const squareType& c
            ) const
            {
                return (a & inv(b)) & c;
            }
    };


private:

    // Private data

        //- Scalar coefficient
        mutable scalarType* scalarCoeffPtr_;

        //- Linear coefficient
        mutable linearType* linearCoeffPtr_;

        //- Square coefficient
        mutable squareType* squareCoeffPtr_;


    // Private Member Functions

        //- Promote to scalar
        scalarType& toScalar();

        //- Promote to linear
        linearType& toLinear();

        //- Promote to square
        squareType& toSquare();


public:

    // Constructors

        //- Construct null
        explicit BlockCoeff();

        //- Construct as copy
        BlockCoeff(const BlockCoeff<Type>&);

        //- Construct from Istream
        BlockCoeff(Istream&);

        //- Clone
        BlockCoeff<Type> clone() const;


    // Destructor

        ~BlockCoeff();

        //- Clear data
        void clear();


    // Member functions

        //- Return active type
        blockCoeffBase::activeLevel activeType() const;

        //- Check pointers: only one type should be active (debug only)
        void checkActive() const;

        // Return as typed.  Fails when asked for the incorrect type

            //- Return as scalar
            const scalarType& asScalar() const;
            scalarType& asScalar();

            //- Return as linear
            const linearType& asLinear() const;
            linearType& asLinear();

            //- Return as square
            const squareType& asSquare() const;
            squareType& asSquare();


        //- Return component
        scalarType component(const direction) const;


    // Member operators

        void operator=(const BlockCoeff<Type>&);


    // IOstream operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const BlockCoeff<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "primitives/BlockCoeff/BlockCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
