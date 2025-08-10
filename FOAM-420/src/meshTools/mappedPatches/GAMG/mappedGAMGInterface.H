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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::mappedGAMGInterface

Description
    Base class for GAMG agglomerated coupled region interface.

SourceFiles
    mappedGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mappedGAMGInterface_H
#define mappedGAMGInterface_H

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "mappedPatches/mappedLduInterface/mappedLduInterface.H"
#include "matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/GAMGAgglomeration/GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class mappedGAMGInterface
:
    public GAMGInterface,
    virtual public mappedLduInterface
{
    // Private data

        //- Reference for the regionCoupledFvPatch from which this is
        //  agglomerated
        const mappedLduInterface& fineMappedLduInterface_;

        //- AMI interface
        autoPtr<AMIPatchToPatchInterpolation> amiPtr_;

        //- Map pointer
        autoPtr<mapDistribute> mapPtr_;

        //- Neighbour GAMG interface index in coarseInterfaces_, or -1 if not found
        mutable label nbrPatchID_;

        //- Done the search for neighbour patch yet?
        mutable bool nbrPatchSearched_;

        //- Size of neighbour coarse patch
        mutable label nbrAddrSize_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mappedGAMGInterface(const mappedGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const mappedGAMGInterface&);

        //- Shift the cell indices so that we get unique values for
        // indices that came from different processors
        void makeIndicesUnique
        (
            labelField& fullRestrictAddressing,
            const labelListList& constructMap
        ) const;

public:

    //- Runtime type information
    TypeName("mappedPatch");

    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        mappedGAMGInterface
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
    virtual ~mappedGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Coupled interface functions

            virtual bool isInterface() const
            {
                return fineMappedLduInterface_.isInterface();
            }

            //- Return neighbour patch number
            virtual label nbrPatchID() const;

            virtual bool owner() const
            {
                return fineMappedLduInterface_.owner();
            }

            virtual const mappedGAMGInterface& nbrPatch() const
            {
                if (nbrPatchID() == -1)
                {
                    FatalErrorInFunction()
                        << "Could not find neighbour patch named" << nl
                        << this->nbrPatchName() << nl << exit(FatalError);
                }
                return dynamic_cast<const mappedGAMGInterface&>
                (
                    coarseInterfaces_[nbrPatchID()]
                );
            }

            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return amiPtr_();
            }

            virtual const mapDistribute& map() const
            {
                return mapPtr_();
            }

            //- Whether we are using AMI or direct map
            virtual bool usingAMI() const
            {
                return fineMappedLduInterface_.usingAMI();
            }

            virtual bool nbrMeshFound() const
            {
                return fineMappedLduInterface_.nbrMeshFound();
            }

            virtual const polyMesh& nbrMesh() const
            {
                NotImplemented;
            }

            virtual const word& patchName() const
            {
                return fineMappedLduInterface_.patchName();
            }

            virtual const word& regionName() const
            {
                return fineMappedLduInterface_.regionName();
            }

            virtual const word& nbrPatchName() const
            {
                return fineMappedLduInterface_.nbrPatchName();
            }

            virtual const word& nbrRegionName() const
            {
                return fineMappedLduInterface_.nbrRegionName();
            }

            template <class Type>
            void interpolate(Field<Type>& pnf) const;

        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                NotImplemented;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatches/GAMG/mappedGAMGInterfaceTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
