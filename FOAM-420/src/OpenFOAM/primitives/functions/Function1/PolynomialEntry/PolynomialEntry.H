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
    (c) 2019-2022 Esi Ltd.

Class
    Foam::Function1Types::PolynomialEntry

Description
    PolynomialEntry container data entry for scalars. Items are stored in a
    list of Tuple2's. Data is input in the form,
    e.g. for an entry \<entryName\> that describes y = x^2 + 2x^3

    \verbatim
        <entryName>  polynomial
        (
            (1    2)
            (2    3)
        );
    \endverbatim

SourceFiles
    PolynomialEntry.C

\*---------------------------------------------------------------------------*/

#ifndef PolynomialEntry_H
#define PolynomialEntry_H

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
                        Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Polynomial
:
    public Function1<Type>
{
    // Private data

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<Type, Type>> coeffs_;

        //- Limits of abscissa value where polynomial is valid
        Pair<scalar> limits_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Polynomial<Type>&);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        Polynomial(const word& entryName, const dictionary& dict);

        //- Construct from components
        Polynomial
        (
            const word& entryName,
            const List<Tuple2<Type, Type>>&,
            const Pair<scalar>& limits=Pair<scalar>(-GREAT, GREAT)
        );

        //- Copy constructor
        Polynomial(const Polynomial& poly);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Polynomial(*this));
        }


    //- Destructor
    virtual ~Polynomial();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return Polynomial value
            Type value(const scalar x) const;

            //- Integrate between two (scalar) values
            Type integrate(const scalar x1, const scalar x2) const;

            //- Integrate function value divided by abscissa between two
            // (scalar) values
            Type integrateYoverX(const scalar x1, const scalar x2) const;

            //- Derivative at the point
            Type derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/PolynomialEntry/PolynomialEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
