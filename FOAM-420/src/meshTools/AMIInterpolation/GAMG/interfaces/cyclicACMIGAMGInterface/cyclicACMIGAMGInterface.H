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
    (c) 2013-2015 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

Class
    Foam::cyclicACMIGAMGInterface

Description
    GAMG agglomerated cyclic ACMI interface.

SourceFiles
    cyclicACMIGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIGAMGInterface_H
#define cyclicACMIGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMILduInterfaceField/cyclicACMILduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicACMIGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicACMILduInterface
{
    // Private data

        //- Reference for the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicACMILduInterface& fineCyclicACMIInterface_;

        //- AMI interface
        autoPtr<AMIPatchToPatchInterpolation> amiPtr_;

        //- Offset from the patch numbers in the original matrix to those in
        // the super-matrix
        label superMatrixPatchOffset_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicACMIGAMGInterface(const cyclicACMIGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicACMIGAMGInterface&);


public:

    //- Runtime type information
    TypeName("cyclicACMI");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicACMIGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );


    //- Destructor
    virtual ~cyclicACMIGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neighbour patch number
            virtual label nbrPatchID() const
            {
                return
                (
                    fineCyclicACMIInterface_.nbrPatchID()
                  + superMatrixPatchOffset_
                );
            }

            virtual bool owner() const
            {
                return fineCyclicACMIInterface_.owner();
            }

            virtual const cyclicACMIGAMGInterface& nbrPatch() const
            {
                return dynamic_cast<const cyclicACMIGAMGInterface&>
                (
                    coarseInterfaces_[nbrPatchID()]
                );
            }

            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return amiPtr_();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return fineCyclicACMIInterface_.transform();
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                //TBD. How to serialise the AMI such that we can stream
                // cyclicACMIGAMGInterface.
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
