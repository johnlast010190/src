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
    Foam::Function1Types::SelectEntry

Description
    Select between different Function1's depending on the value of the input.
    A list of N values is supplied and a list of N+1 dictionaries containing
    Function1s to be used for each range of values.

    \verbatim
        <entryName>  select ( 300 600 )
        (
            {
                (first function - used for values < 300)
            }
            {
                (second function - used for values > 300 and < 600)
            }
            {
                (third function - used for values > 600)
            }
        );
    \endverbatim

SourceFiles
    SelectEntry.C

\*---------------------------------------------------------------------------*/

#ifndef SelectEntry_H
#define SelectEntry_H

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
                        Class Select Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Select
:
    public Function1<Type>
{
    // Private data

        //- List of values at which we switch between the child Function1s
        List<scalar> switchValues_;

        //- Child functions, one more than the number of switch values
        PtrList<Function1<Type>> childFunctions_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Select<Type>&);

        void checkValid() const;

public:

    //- Runtime type information
    TypeName("select");


    // Constructors

        Select(const word& entryName, const dictionary& dict);

        //- Construct from components
        Select
        (
            const word& entryName,
            const List<scalar>& switchValues,
            const PtrList<Function1<Type>>& childFunctions
        );

        //- Copy constructor
        Select(const Select& sel);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Select(*this));
        }


    //- Destructor
    virtual ~Select();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return value
            Type value(const scalar x) const;

            //- Integrate between two (scalar) values
            Type integrate(const scalar x1, const scalar x2) const;

            //- Integrate function value divided by abscissa between two
            // (scalar) values
            Type integrateYoverX(const scalar x1, const scalar x2) const;

            //- Derivative
            virtual Type derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/SelectEntry/SelectEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
