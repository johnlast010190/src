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
    Foam::Function1Types::powerLaw

Description
    Scalar power law function with support for an offset level.

        \f[
            a (t - b)^c + d
        \f]

    where

    \vartable
        symbol               | Description       | Data type
        a                    | Multiplicator     | scalar
        b                    | Value offset      | scalar
        c                    | Exponent          | scalar
        d                    | Variable          | scalar
        maxValue             | Optional maximum
        forcePositiveBracket | Recompute minimum x so that output is non-negative
    \endvartable

    Example for a scalar:
    \verbatim
        <entryName> powerLaw;
        <entryName>Coeffs
        {
            a        10;
            b        273;
            c        0.1;
            d        10;
            maxValue 10;
            forcePositiveBracket true;
        }
    \endverbatim

SourceFiles
    powerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public Function1<scalar>
{
    // Private data

        //- Coefficient 1
        scalar a_;

        //- Coefficient 2
        scalar b_;

        //- Coefficient 3
        scalar c_;

        //- Coefficient 4
        scalar d_;

        //- Maximum value
        scalar maxValue_;

        //- Minimum value
        scalar minValue_;

        //- Clip?
        Switch clip_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);

        //- Disallow default bitwise assignment
        void operator=(const powerLaw&);


public:

    // Runtime type information
    TypeName("powerLaw");


    // Constructors

        //- Construct from entry name and dictionary
        powerLaw
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        powerLaw(const powerLaw& se);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new powerLaw(*this));
        }


    //- Destructor
    virtual ~powerLaw();


    // Member Functions

        //- Return value for input variable
        scalar value(const scalar x) const;

        //- Derivative
        scalar derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
