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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalCyclicPolyPatch

Description
    Non-conformal cyclic poly patch. As nonConformalCoupledPolyPatch, but the
    neighbouring patch is local and known and is made available by this class.

See also
    Foam::nonConformalCoupledPolyPatch

SourceFiles
    nonConformalCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicPolyPatch_H
#define nonConformalCyclicPolyPatch_H

#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.H"
#include "patchToPatch/intersection/intersectionPatchToPatch.H"
#include "patchToPatch/rays/raysPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class nonConformalCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCyclicPolyPatch
:
    public cyclicPolyPatch,
    public nonConformalCoupledPolyPatch
{
protected:

    // Protected data

        //- Is the intersection engine up to date?
        mutable bool intersectionIsValid_;

        //- Patch-to-patch intersection engine
        mutable patchToPatches::intersection intersection_;

        //- Is the intersection engine up to date?
        mutable bool raysIsValid_;

        //- Patch-to-patch intersection engine
        mutable patchToPatches::rays rays_;


    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) override;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) override;

        //- Initialise the patches for moving points
        virtual void initMovePoints
        (
            PstreamBuffers& pBufs,
            const pointField&
        ) override;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&) override;

        //- Clear geometry
        virtual void clearGeom() override;

        //- Reset the patch name
        virtual void rename(const wordList& newNames) override;

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex) override;


public:

    //- Runtime type information
    TypeName("nonConformalCyclic");


    // Constructors

        //- Construct from components
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& nbrPatchName,
            const word& origPatchName,
            const cyclicTransform& transform=cyclicTransform(true)
        );

        //- Construct from dictionary
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        nonConformalCyclicPolyPatch
        (
            const nonConformalCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalCyclicPolyPatch
        (
            const nonConformalCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName,
            const word& origPatchName
        );

        //- Construct given the original patch and a map
        nonConformalCyclicPolyPatch
        (
            const nonConformalCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm
        ) const override
        {
            return autoPtr<polyPatch>
            (
                new nonConformalCyclicPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const override
        {
            return autoPtr<directPolyPatch>
            (
                new nonConformalCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName(),
                    origPatchName()
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const override
        {
            return autoPtr<directPolyPatch>
            (
                new nonConformalCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~nonConformalCyclicPolyPatch();


    // Member Functions

        //- Neighbour patch
        const nonConformalCyclicPolyPatch& nbrPatch() const
        {
            return
                refCast<const nonConformalCyclicPolyPatch>
                (
                    cyclicPolyPatch::nbrPatch()
                );
        }

        //- Inherit the cyclic owner method
        using cyclicPolyPatch::owner;

        //- Inherit the cyclic neighbour method
        using cyclicPolyPatch::neighbour;

        //- Inherit the cyclic transform method
        using cyclicPolyPatch::transform;

        //- Is this patch coupled? Returns false. For NCC patches, the poly
        //  mesh is considered non-coupled whilst the finite volume mesh is
        //  considered coupled.
        virtual bool coupled() const override
        {
            return false;
        }

        //- Access the intersection engine
        const patchToPatches::intersection& intersection() const;

        //- Access the rays engine
        const patchToPatches::rays& rays() const;

        //- Compute a ray intersection across the coupling
        patchToPatch::procFace ray
        (
            const scalar fraction,
            const label origFacei,
            const vector& p,
            const vector& n,
            point& nbrP
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
