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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2015 Esi Ltd.

Class
    Foam::offsetGAMGInterface

Description
    GAMG agglomerated interface.

SourceFiles
    offsetGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef offsetGAMGInterface_H
#define offsetGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "offsetLduInterface/offsetLduInterface.H"
#include "offsetLduInterface/offsetLduInterfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class offsetGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class offsetGAMGInterface
:
    public GAMGInterface,
    virtual public offsetLduInterfaceBase
{
    // Private data members

        autoPtr< GAMGInterface > slaveGAMGInterface_;

        mutable lduInterfacePtrsList slaveCoarseInterfaces_;

    // Private Member Functions

        //- Coarse interfaces are not populated at construction, so this
        // function is always called to make sure they have been obtained
        void grabSlaveCoarseInterfaces() const;


    //- Disallow default bitwise copy construct
        offsetGAMGInterface(const offsetGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const offsetGAMGInterface&);


public:

    //- Runtime type information
    TypeName("offset");

    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        offsetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from Istream
        offsetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~offsetGAMGInterface();


    // Member Functions

        // Access

            //- Return underlying interface we are controlling
            virtual const lduInterface& slaveInterface() const
            {
                grabSlaveCoarseInterfaces();
                return slaveGAMGInterface_();
            }

            const GAMGInterface& slaveGAMGInterface() const
            {
                grabSlaveCoarseInterfaces();
                return slaveGAMGInterface_();
            }

            GAMGInterface& slaveGAMGInterface()
            {
                grabSlaveCoarseInterfaces();
                return slaveGAMGInterface_();
            }

    //- Return size
            virtual label size() const
            {
                return slaveGAMGInterface().size();
            }

            virtual label index() const
            {
                return slaveGAMGInterface().index();
            }

            virtual const lduInterfacePtrsList& coarseInterfaces() const
            {
                return slaveGAMGInterface().coarseInterfaces();
            }

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return slaveGAMGInterface().faceCells();
            }

            //- Return (local)face restrict addressing
            virtual const labelList& faceRestrictAddressing() const
            {
                return slaveGAMGInterface().faceRestrictAddressing();
            }

            //- Return non-const access to face restrict addressing
            virtual labelList& faceRestrictAddressing()
            {
                return slaveGAMGInterface().faceRestrictAddressing();
            }

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;


        // Interface transfer functions

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        // Agglomeration

            //- Agglomerate the given fine-level coefficients and return
            virtual tmp<scalarField> agglomerateCoeffs
            (
                const scalarField& fineCoeffs
            ) const;


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
