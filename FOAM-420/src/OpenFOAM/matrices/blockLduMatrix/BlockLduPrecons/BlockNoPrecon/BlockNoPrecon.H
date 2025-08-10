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
    BlockNoPrecon

Description
    No preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockNoPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockNoPrecon_H
#define BlockNoPrecon_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class BlockLduPrecon;

/*---------------------------------------------------------------------------*\
                       Class BlockNoPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockNoPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockNoPrecon(const BlockNoPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockNoPrecon&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        BlockNoPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const word fieldName
        )
        :
            BlockLduPrecon<Type>(matrix)
        {}


    // Destructor

        virtual ~BlockNoPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const
        {
            x = b;
        }

        //- Execute preconditioning
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const
        {
            precondition(xT, bT);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
