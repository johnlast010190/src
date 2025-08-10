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
    (c) 2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::Function2s::None

Description
    Undefined Function2 which returns an error when called.

\*---------------------------------------------------------------------------*/

#ifndef None2_H
#define None2_H

#include "primitives/functions/Function2/Function2/Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class None Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class None
:
    public FieldFunction2<Type, None<Type>>
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        fileName dictName_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from name and dictionary
        None(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~None();


    // Member Functions

        //- Generates an error
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const
        {
            return tmp<Function2<Type>>(new None(*this));
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const None<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2
} // End namespace Foam

#ifdef NoRepository
    #include "primitives/functions/Function2/None/None2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
