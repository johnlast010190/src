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
    (c) 2016-2017 Wikki Ltd
    (c) 2018-2022 Esi Ltd.

Class
    Foam::cyclicFaPatch

Description
    Cyclic-plane patch.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaPatch_H
#define cyclicFaPatch_H

#include "faMesh/faPatches/basic/coupled/coupledFaPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/cyclicLduInterface.H"
#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFaPatch
:
    public coupledFaPatch,
    public cyclicLduInterface
{
    // Private data
    scalar matchTolerance_;

    // Private member functions

    void calcTransforms();

protected:

    // Protected static data

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from dictionary
        cyclicFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm),
            matchTolerance_(1e-4)
        {}


    // Destructor

        virtual ~cyclicFaPatch()
        {}


    // Member functions

        // Access

            //- Return neighbour
               virtual label nbrPatchID() const
               {
                   return coupledFaPatch::ngbPolyPatchIndex();
               }

               virtual bool owner() const
               {
                   return 1;
               }

               //- Return processor number
               virtual const cyclicLduInterface& nbrPatch() const
               {
                   return refCast<const cyclicFaPatch>
                   (
                       this->boundaryMesh()
                   );
               }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return coupledFaPatch::forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return coupledFaPatch::reverseT();
            }

            //- Return dummy null transform
            virtual const transformer& transform() const
            {
                return transformer::null;
            }

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry(PstreamBuffers&);

            //- Calculate the patch geometry
            virtual void calcGeometry(PstreamBuffers&);

            //- Initialise the patches for moving points
            virtual void initMovePoints(PstreamBuffers&, const pointField&);

            //- Correct patches after moving points
            virtual void movePoints(PstreamBuffers&, const pointField&);

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
