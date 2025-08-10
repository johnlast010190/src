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
    (c) 2020-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::Function2s::OneConstant

Description
    Templated function of two variables that returns the corresponding 1 (one).

Usage
    \verbatim
        <name> one;
    \endverbatim

SourceFiles
    OneConstant.C

\*---------------------------------------------------------------------------*/

#ifndef OneConstant2_H
#define OneConstant2_H

#include "primitives/functions/Function2/Function2/Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                         Class OneConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class OneConstant
:
    public FieldFunction2<Type, OneConstant<Type>>
{

public:

    // Runtime type information
    TypeName("one");


    // Constructors

        //- Construct from name
        OneConstant(const word& name);

        //- Construct from name and dictionary
        OneConstant(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~OneConstant();


    // Member Functions

        //- Return One
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Return One field
        virtual inline tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const
        {
            return tmp<Function2<Type>>(new OneConstant(*this));
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OneConstant<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OneConstant2I.H"

#ifdef NoRepository
    #include "primitives/functions/Function2/One/OneConstant2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
