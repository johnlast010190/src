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
    (c) 2022 Esi Ltd.

Class
    Foam::processorLduInterface

Description
    An abstract base class for processor coupled interfaces.

SourceFiles
    processorLduInterface.C
    processorLduInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorLduInterface_H
#define processorLduInterface_H

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "primitives/transform/transformer/transformer.H"
#include "fields/Fields/primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorLduInterface Declaration
\*---------------------------------------------------------------------------*/

class processorLduInterface
{
    // Private data

        //- Send buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> sendBuf_;

        //- Receive buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> receiveBuf_;

        //- Resize the buffer if required
        void resizeBuf(List<char>& buf, const label size) const;


public:

    //- Runtime type information
    TypeName("processorLduInterface");


    // Constructors

        //- Construct null
        processorLduInterface();


    //- Destructor
    virtual ~processorLduInterface();


    // Member Functions

        // Access

            //- Return communicator used for parallel communication
            virtual label comm() const = 0;

            //- Return processor number (rank in communicator)
            virtual int myProcNo() const = 0;

            //- Return neigbour processor number (rank in communicator)
            virtual int neighbProcNo() const = 0;

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return message tag used for sending
            virtual int tag() const = 0;


        // Transfer functions

            //- Raw send function
            template<class Type>
            void send
            (
                const Pstream::commsTypes commsType,
                const UList<Type>&
            ) const;

            //- Raw field receive function
            template<class Type>
            void receive
            (
                const Pstream::commsTypes commsType,
                UList<Type>&
            ) const;

            //- Raw field receive function returning field
            template<class Type>
            tmp<Field<Type>> receive
            (
                const Pstream::commsTypes commsType,
                const label size
            ) const;


            //- Raw field send function with data compression
            template<class Type>
            void compressedSend
            (
                const Pstream::commsTypes commsType,
                const UList<Type>&
            ) const;

            //- Raw field receive function with data compression
            template<class Type>
            void compressedReceive
            (
                const Pstream::commsTypes commsType,
                UList<Type>&
            ) const;

            //- Raw field receive function with data compression returning field
            template<class Type>
            tmp<Field<Type>> compressedReceive
            (
                const Pstream::commsTypes commsType,
                const label size
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/lduMatrix/lduAddressing/lduInterface/processorLduInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
