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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::cyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    cyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFvPatch_H
#define cyclicFvPatch_H

#include "fvMesh/fvPatches/basic/coupled/coupledFvPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/cyclicLduInterface.H"
#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFvPatch
:
    public coupledFvPatch,
    public cyclicLduInterface
{
    // Private data

        const cyclicPolyPatch& cyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            cyclicPolyPatch_(refCast<const cyclicPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return true if this patch field is constrained
            virtual bool constrained() const
            {
                return true;
            }

            //- Return local reference cast into the cyclic patch
            const cyclicPolyPatch& cyclicPatch() const
            {
                return cyclicPolyPatch_;
            }

            //- Return neighbour
            virtual label nbrPatchID() const
            {
                return cyclicPolyPatch_.nbrPatchID();
            }

            virtual bool owner() const
            {
                return cyclicPolyPatch_.owner();
            }

            //- Return processor number
            virtual const cyclicFvPatch& nbrPatch() const
            {
                return refCast<const cyclicFvPatch>
                (
                    this->boundaryMesh()[cyclicPolyPatch_.nbrPatchID()]
                );
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicPolyPatch_.transform();
            }

            const cyclicFvPatch& neighbFvPatch() const
            {
                return refCast<const cyclicFvPatch>
                (
                    this->boundaryMesh()[cyclicPolyPatch_.nbrPatchID()]
                );
            }


            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
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
