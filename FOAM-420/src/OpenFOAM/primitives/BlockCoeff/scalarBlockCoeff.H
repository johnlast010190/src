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

Class
    BlockCoeff specialisation for terminal class.

Description
    Specialisation of a block coefficient for a scalar: always a scalae.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockCoeff_H
#define scalarBlockCoeff_H

#include "primitives/BlockCoeff/BlockCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

template<>
class BlockCoeff<scalar>;
Ostream& operator<<(Ostream&, const BlockCoeff<scalar>&);

/*---------------------------------------------------------------------------*\
                          Class BlockCoeff Declaration
\*---------------------------------------------------------------------------*/

template<>
class BlockCoeff<scalar>
:
    public blockCoeffBase
{
public:

    // Public data types

        //- Component type
        typedef scalar Type;

        typedef Type xType;

        //- Coefficient type
        typedef scalar scalarType;


    //- Multiplication trait
    class multiply
    {
    public:

        multiply() {}

        Type operator()(const scalarType& c, const Type& x) const
        {
            return c*x;
        }
        Type operator()(const BlockCoeff<Type>& c, const Type& x) const
        {
            return c.asScalar()*x;
        }

        // Inverse functions
        scalarType inverse(const scalarType& c) const
        {
            return 1.0/c;
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
    };


private:

    // Private data

        //- Scalar coefficient
        scalar scalarCoeff_;


public:

    // Constructors

        //- Construct null
        explicit BlockCoeff();

        //- Construct as copy
        BlockCoeff(const BlockCoeff<scalar>&);

        //- Construct from Istream
        BlockCoeff(Istream&);

        //- Clone
        BlockCoeff<scalar> clone() const;


    // Destructor

        ~BlockCoeff();


    // Member functions

        //- Return as scalar
        const scalar& asScalar() const
        {
            return scalarCoeff_;
        }

        scalar& asScalar()
        {
            return scalarCoeff_;
        }


        //- Return component
        scalar component(const direction) const;


    // Member operators

        void operator=(const BlockCoeff<scalar>&);


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const BlockCoeff<scalar>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
