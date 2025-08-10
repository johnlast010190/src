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
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicGAMGInterface

Description
    GAMG agglomerated cyclic interface.

SourceFiles
    cyclicGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGAMGInterface_H
#define cyclicGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/cyclicLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Neigbour patch number
        label nbrPatchID_;

        //- Am I owner?
        bool owner_;

        //- Transformation
        transformer transform_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicGAMGInterface(const cyclicGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicGAMGInterface&);


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicGAMGInterface
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
        cyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~cyclicGAMGInterface();


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
                return nbrPatchID_;
            }

            virtual bool owner() const
            {
                return owner_;
            }

            virtual const cyclicGAMGInterface& nbrPatch() const
            {
                return dynamic_cast<const cyclicGAMGInterface&>
                (
                    coarseInterfaces_[nbrPatchID()]
                );
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return transform_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
