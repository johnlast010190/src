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
    (c) 2017 OpenFOAM Foundation
    (c) 2020-2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::Function1Types::ZeroConstant

Description
    Templated function that returns the corresponding 0 (zero).

    Usage:
    \verbatim
        <entryName> zero;
    \endverbatim

SourceFiles
    ZeroConstant.C

\*---------------------------------------------------------------------------*/

#ifndef ZeroConstant_H
#define ZeroConstant_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class ZeroConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ZeroConstant
:
    public Function1<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ZeroConstant<Type>&) = delete;


public:

    // Runtime type information
    TypeName("zero");


    // Constructors

        //- Construct from entry name and dictionary
        ZeroConstant(const word& entryName, const dictionary& dict);

        //- Construct and return a clzero
        virtual tmp<Function1<Type>> clzero() const
        {
            return tmp<Function1<Type>>(new ZeroConstant<Type>(*this));
        }


    //- Destructor
    virtual ~ZeroConstant();


    // Member Functions

        //- Return constant value
        Type value(const scalar) const;

        //- Integrate between two values
        Type integrate(const scalar x1, const scalar x2) const;

        //- Derivative
        Type derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/Zero/ZeroConstant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
