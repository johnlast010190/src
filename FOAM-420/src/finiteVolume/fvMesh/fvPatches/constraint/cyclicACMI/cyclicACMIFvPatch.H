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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicACMIFvPatch

Description
    Cyclic patch for Arbitrarily Coupled Mesh Interface (ACMI)

SourceFiles
    cyclicACMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIFvPatch_H
#define cyclicACMIFvPatch_H

#include "fvMesh/fvPatches/basic/coupled/coupledFvPatch.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMILduInterfaceField/cyclicACMILduInterface.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMIPolyPatch/cyclicACMIPolyPatch.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicACMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIFvPatch
:
    public coupledFvPatch,
    public cyclicACMILduInterface
{
    // Private data

        const cyclicACMIPolyPatch& cyclicACMIPolyPatch_;


protected:

    // Protected Member functions

        //- Update the patch areas after AMI update
        void updateAreas() const;

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicACMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            cyclicACMILduInterface(),
            cyclicACMIPolyPatch_(refCast<const cyclicACMIPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicACMIPolyPatch& cyclicACMIPatch() const
            {
                return cyclicACMIPolyPatch_;
            }

            //- Return neighbour
            virtual label nbrPatchID() const
            {
                return cyclicACMIPolyPatch_.nbrPatchID();
            }

            virtual bool owner() const
            {
                return cyclicACMIPolyPatch_.owner();
            }

            //- Return neighbour fvPatch
            virtual const cyclicACMIFvPatch& nbrPatch() const
            {
                return refCast<const cyclicACMIFvPatch>
                (
                    this->boundaryMesh()[cyclicACMIPolyPatch_.nbrPatchID()]
                );
            }

            //- Return neighbour
            virtual label nonOverlapPatchID() const
            {
                return cyclicACMIPolyPatch_.nonOverlapPatchID();
            }

            //- Return non-overlapping fvPatch
            virtual const fvPatch& nonOverlapPatch() const
            {
                return this->boundaryMesh()[nonOverlapPatchID()];
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                const AMIPatchToPatchInterpolation& AMI =
                    cyclicACMIPolyPatch_.AMI();

                updateAreas();

                return AMI;
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicACMIPolyPatch_.transform();
            }

            const cyclicACMIFvPatch& neighbFvPatch() const
            {
                return refCast<const cyclicACMIFvPatch>
                (
                    this->boundaryMesh()[cyclicACMIPolyPatch_.nbrPatchID()]
                );
            }

            //- Return true if this patch is coupled. This is equivalent
            //  to the coupledPolyPatch::coupled() if parallel running or
            //  both sides present, false otherwise
            virtual bool coupled() const;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            //- Interpolate (make sure to have uptodate areas)
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld
            ) const
            {
                updateAreas();

                return
                    cyclicACMIPolyPatch_.cyclicAMIPolyPatch::interpolate
                    (
                        fld
                    );
            }

            //- Interpolate (make sure to have uptodate areas)
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tfld
            ) const
            {
                return interpolate(tfld());
            }


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
