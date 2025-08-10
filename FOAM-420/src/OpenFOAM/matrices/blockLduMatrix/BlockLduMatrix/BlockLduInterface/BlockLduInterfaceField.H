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
    (c) 2004-2006 Ivor Clifford

Class
    BlockLduInterfaceField

Description
    Abstract base class for interface fields with block coefficients

Author
    Ivor Clifford, ivor.clifford@gmail.com

SourceFiles
    BlockLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterfaceField_H
#define BlockLduInterfaceField_H

#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
class CoeffField;

/*---------------------------------------------------------------------------*\
                           Class BlockLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterfaceField
:
    public lduInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduInterfaceField(const BlockLduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterfaceField&);


public:

    //- Runtime Type information
    TypeName("BlockLduInterfaceField");

    // Constructors

        //- Construct given coupled patch
        BlockLduInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}


    // Destructor

        virtual ~BlockLduInterfaceField()
        {}



    // Member Functions

        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
