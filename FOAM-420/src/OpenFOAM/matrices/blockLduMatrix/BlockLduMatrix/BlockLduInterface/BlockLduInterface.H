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
    (c) 2004-2011 H. Jasak

Class
    BlockLduInterface

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterface_H
#define BlockLduInterface_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/FieldFields/FieldField/FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

/*---------------------------------------------------------------------------*\
                     Class BlockLduInterface Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterface
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduInterface(const BlockLduInterface&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterface&);


public:

    //- Runtime type information
    TypeName("BlockLduInterface");


    // Destructor

        virtual ~BlockLduInterface()
        {}


    // Member Functions

        // Coupled interface matrix update

            //- Initialise matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const BlockLduMatrix<Type>& matrix,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const
            {}

            //- Update result based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const BlockLduMatrix<Type>& matrix,
                const FieldField<CoeffField, Type>& coeffs,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
