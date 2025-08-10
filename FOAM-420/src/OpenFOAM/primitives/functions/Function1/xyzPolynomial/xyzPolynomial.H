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
    (c) 2022 Esi Ltd.

Class
    Foam::Function1Types::xyzPolynomial

Description
    xyzPolynomial function that takes list of values and exponents and creates
    generic three variable variable polynomial.

    Example:
    <entryName> = 300 + 0.5x + 3y + 0.3xy^0.5

    \verbatim
        <entryName> xyzPolynomial
        (
            (300 (0 0 0))   // 300
            (0.5 (1 0 0))   // 0.5*x
            (3   (0 1 0))   // 3*y
            (0.3 (1 0.5 0)) // 0.3*x*y^0.5
        );

    \endverbatim

SourceFiles
    xyzPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef xyzPolynomial_H
#define xyzPolynomial_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/Pair/Pair.H"
#include "primitives/Tuple2/Tuple2.H"
#include "primitives/functions/Function1/Function1/Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class xyzPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class xyzPolynomial
:
    public Function1<Type>
{
    // Private data

        //- Component vector field
        List<Tuple2<Type, vector>> inputComp_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const xyzPolynomial<Type>&);


public:

    //- Runtime type information
    TypeName("xyzPolynomial");


    // Constructors

        xyzPolynomial(const word& entryName, const dictionary& dict);

        //- Construct from components
        xyzPolynomial
        (
            const word& entryName,
            const List<Tuple2<Type, vector>>& coeffs
        );

        //- Copy constructor
        xyzPolynomial(const xyzPolynomial& poly);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const override
        {
            return tmp<Function1<Type>>(new xyzPolynomial(*this));
        }


    //- Destructor
    virtual ~xyzPolynomial();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t) override;


        // Evaluation

            //- Return xyzPolynomial value
            Type value(const vector& xyz) const override;

            //- Integrate between two (scalar) values
            Type integrate(const scalar x1, const scalar x2) const override;

            //- Integrate function value divided by abscissa between two
            // (scalar) values
            Type integrateYoverX(const scalar x1, const scalar x2) const
            override;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const override;

        //- Read function 1
        void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/xyzPolynomial/xyzPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
