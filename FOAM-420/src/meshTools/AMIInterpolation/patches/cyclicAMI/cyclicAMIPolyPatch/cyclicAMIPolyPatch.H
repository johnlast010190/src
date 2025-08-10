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
    (c) 2020 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

Class
    Foam::cyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPolyPatch_H
#define cyclicAMIPolyPatch_H

#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicTransform.H"
#include "AMIInterpolation/AMIInterpolation/AMIPatchToPatchInterpolation.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/Identifiers/patch/coupleGroupIdentifier.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPolyPatch
:
    public coupledPolyPatch,
    public cyclicTransform
{
protected:

    // Protected data

        //- Name of cyclic neighbour patch
        mutable word nbrPatchName_;

        //- Optional patchGroup to find neighbour patch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of cyclic neighbour patch
        mutable label nbrPatchID_;

        //- AMI interpolation class
        mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;

        //- AMI method
        const AMIPatchToPatchInterpolation::interpolationMethod AMIMethod_;

        //- Flag to indicate that slave patch should be reversed for AMI
        const bool AMIReverse_;

        //- Flag to indicate that patches should match/overlap
        bool AMIRequireMatch_;

        //- Low weight correction threshold for AMI
        scalar AMILowWeightCorrection_;

        //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;

        // Change of topology as AMI is updated

            //- Flag to indicate that new AMI faces will created
            //  Set by the call to changeTopology
            mutable bool createAMIFaces_;

            //- Move face centres (default = no)
            bool moveFaceCentres_;

            labelListList srcFaceIDs_;

            labelListList tgtFaceIDs_;

            //- Temporary storage for AMI face areas
            mutable vectorField faceAreas0_;

            //- Temporary storage for AMI face centres
            mutable vectorField faceCentres0_;


    // Protected Member Functions

        // Topology change

            //- Collect faces to remove in the topoChange container
            virtual bool removeAMIFaces(polyTopoChange& topoChange);

            //- Collect faces to add in the topoChange container
            virtual bool addAMIFaces(polyTopoChange& topoChange);

            //- Set properties of newly inserted faces after topological changes
            virtual void setAMIFaces();

            //- Helper to re-apply the geometric scaling lost during mesh
            //- updates
            virtual void restoreScaledGeometry();

        //- Reset the AMI interpolator, supply patch points
//        virtual void resetAMI(const UList<point>& points) const;

        //- Reset the AMI interpolator

        virtual void resetAMI
        (
            const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
        ) const;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

        //- Clear geometry
        virtual void clearGeom();

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex);


public:

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        cyclicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicAMIPolyPatch(const cyclicAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new cyclicAMIPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new cyclicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
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
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new cyclicAMIPolyPatch
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
    virtual ~cyclicAMIPolyPatch();


    // Member Functions

        // Access

            //- Tolerance used e.g. for area calculations/limits
            static const scalar tolerance_;

            //- Flag to indicate whether the AMI can be reset
            inline bool canResetAMI() const;

            //- Return access to the createAMIFaces flag
            inline bool createAMIFaces() const;

            //- Return true if this patch changes the mesh topology
            //  True when createAMIFaces is true
            virtual bool changeTopology() const;

            //- Set topology changes in the polyTopoChange object
            virtual bool setTopology(polyTopoChange& topoChange);

            //- Is patch 'coupled'. Note that on AMI the geometry is not
            //  coupled but the fields are!
            virtual bool coupled() const
            {
                return false;
            }

            virtual bool implicitNonCovered() const
            {
                return true;
            }

            //- Neighbour patch name
            inline const word& nbrPatchName() const;

            //- Neighbour patch ID
            virtual label nbrPatchID() const;

            //- Does this side own the patch?
            virtual bool owner() const;

            //- Return a reference to the neighbour patch
            virtual const cyclicAMIPolyPatch& nbrPatch() const;

            //- Return a reference to the projection surface
            const autoPtr<searchableSurface>& surfPtr() const;

            //- Return a reference to the AMI interpolator
            const AMIPatchToPatchInterpolation& AMI() const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;

            //- total weight access
            const scalarField& lowCorrAMIWeights() const;

            //- Return access to the initial face areas
            //  Used for topology change
            inline vectorField& faceAreas0() const;

            //- Return access to the initial face centres
            //  Used for topology change
            inline vectorField& faceCentres0() const;


        // Transformations

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicTransform::transform();
            }


        // Interpolations

            //- Interpolate field
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const;

            //- Interpolate tmp field
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const;

            //- Low-level interpolate List
            template<class Type, class CombineOp>
            void interpolate
            (
                const UList<Type>& fld,
                const CombineOp& cop,
                List<Type>& result,
                const UList<Type>& defaultValues = UList<Type>()
            ) const;


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Return face index on neighbour patch which shares point p
        //  following trajectory vector n
        label pointFace
        (
            const label facei,
            const vector& n,
            point& p
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
