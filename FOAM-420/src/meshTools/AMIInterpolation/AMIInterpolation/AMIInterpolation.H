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
    (c) 2016-2020 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

Class
    Foam::AMIInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitive patches with an arbitrary mesh interface (AMI).

    Based on the algorithm given in:

        Conservative interpolation between volume meshes by local Galerkin
        projection, Farrell PE and Maddison JR, 2011, Comput. Methods Appl.
        Mech Engrg, Volume 200, Issues 1-4, pp 89-100

    Interpolation requires that the two patches should have opposite
    orientations (opposite normals).  The 'reverseTarget' flag can be used to
    reverse the orientation of the target patch.


SourceFiles
    AMIInterpolation.C
    AMIInterpolationName.C
    AMIInterpolationParallelOps.C

\*---------------------------------------------------------------------------*/

#ifndef AMIInterpolation_H
#define AMIInterpolation_H

#include "db/typeInfo/className.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "meshes/treeBoundBox/treeBoundBoxList.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "AMIInterpolation/faceAreaIntersect/faceAreaIntersect.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "primitives/ops/ops.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/primitiveShapes/point/pointList.H"
#include "primitives/transform/transformer/transformer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AMIInterpolationName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(AMIInterpolation);


/*---------------------------------------------------------------------------*\
                      Class AMIInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class AMIInterpolation
:
    public AMIInterpolationName
{
public:

    // Public data
        static const scalar relTolDef_;
        static const scalar tolDef_;
        static const scalar cosMatchAngleDef_;
        static const scalar maxAMIWeightScaleDef_;
        static const bool outputDef_;

        //- Enumeration specifying interpolation method
        enum interpolationMethod
        {
            imDirect,
            imMapNearest,
            imFaceAreaWeight,
            imPartialFaceAreaWeight
        };

        //- Convert interpolationMethod to word representation
        static word interpolationMethodToWord
        (
            const interpolationMethod& method
        );

        //- Convert word to interpolationMethod
        static interpolationMethod wordTointerpolationMethod
        (
            const word& method
        );

        //- Calculate the patch face magnitudes for the given tri-mode
        template
        <
            class pType
        >
        static tmp<scalarField> patchMagSf
        (
            const pType& patch,
            const faceAreaIntersect::triangulationMode triMode
        );


private:

    // Private data

        //- Interpolation method
        const word methodName_;

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Flag to indicate that the two patches must be matched/an overlap
        //  exists between them
        const bool requireMatch_;

        //- Index of processor that holds all of both sides. -1 in all other
        //  cases
        label singlePatchProc_;

        //- Threshold weight below which interpolation is deactivated
        scalar lowWeightCorrection_;

        //- Cosine of minimum angle required for face match
        const scalar cosMatchAngle_;

        //- Max allowed scale of the weights
        const scalar maxScale_;

        //- Relative matching distance tolerance as a fraction of bounding box
        const scalar relTol_;

        //- Absolute matching distance tolerance (i.e. max possible distance
        // between matched faces in m)
        const scalar tol_;

        //- Flag to control reporting
        const bool output_;

        // Source patch

            //- Source face areas
            scalarField srcMagSf_;

            //- Source face normals
            vectorField srcNf_;

            //- Addresses of target faces per source face
            labelListList srcAddress_;

            //- Weights of target faces per source face
            scalarListList srcWeights_;

            //- Weights of target faces per source face
            scalarList srcSymWeights_;

            //- Sum of weights of target faces per source face
            scalarField srcWeightsSum_;

            //- Centroid of target faces per source face
            pointListList srcCentroids_;

            //- Transformation map per source face
            labelListList srcTransforms_;

        // Target patch

            //- Target face areas
            scalarField tgtMagSf_;

            //- Target face normals
            vectorField tgtNf_;

            //- Addresses of source faces per target face
            labelListList tgtAddress_;

            //- Weights of source faces per target face
            scalarListList tgtWeights_;

            //- Weights of source faces per target face
            scalarList tgtSymWeights_;

            //- Sum of weights of source faces per target face
            scalarField tgtWeightsSum_;

            //- Centroid of source faces per target face
            pointListList tgtCentroids_;

            //- Transformation map per target face
            labelListList tgtTransforms_;

            //- Source and target addressings for overlap and nonOverlap faces
            mutable labelList* overlapSrcAddrPtr_;

            mutable labelList* nonOverlapSrcAddrPtr_;

            mutable labelList* overlapTgtAddrPtr_;

            mutable labelList* nonOverlapTgtAddrPtr_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;

        //- Source map pointer - parallel running only
        autoPtr<mapDistribute> srcMapPtr_;

        //- Target map pointer - parallel running only
        autoPtr<mapDistribute> tgtMapPtr_;

        //- Source transformation per rotation (used by cyclicPeriodicAMI)
        List<transformer> srcAMITransforms_;

        //- Target transformation per rotation (used by cyclicPeriodicAMI)
        List<transformer> tgtAMITransforms_;

        //- Up-to-date flag
        bool upToDate_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        AMIInterpolation(const AMIInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const AMIInterpolation&);


        // Parallel functionality

            //- Calculate if patches are on multiple processors
            label calcDistribution
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            ) const;

            label calcOverlappingProcs
            (
                const List<treeBoundBoxList>& procBb,
                const treeBoundBox& bb,
                boolList& overlaps
            ) const;

            void distributePatches
            (
                const mapDistribute& map,
                const TargetPatch& pp,
                const globalIndex& gi,
                List<faceList>& faces,
                List<pointField>& points,
                List<labelList>& tgtFaceIDs
            ) const;

            void distributeAndMergePatches
            (
                const mapDistribute& map,
                const TargetPatch& tgtPatch,
                const globalIndex& gi,
                faceList& tgtFaces,
                pointField& tgtPoints,
                labelList& tgtFaceIDs
            ) const;

            autoPtr<mapDistribute> calcProcMap
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            ) const;


        // Initialisation

            //- Project points to surface
            void projectPointsToSurface
            (
                const searchableSurface& surf,
                pointField& pts
            ) const;


        // Debug weights

            //- creates VTK with the source and target patches
            //  In each surface the sumWeights and symWeights overlapping areas
            //  normalised with the face areas are written
            //  NOTE: the symWeights are not modified in the normalisation
            //  function. So they are the same before and after


        // Evaluation

            //- Normalise the (area) weights - suppresses numerical error in
            //  weights calculation
            //  NOTE: if area weights are incorrect by 'a significant amount'
            //     normalisation may stabilise the solution, but will introduce
            //     numerical error!
            static void normaliseWeights
            (
                const scalarField& patchAreas,
                const word& patchName,
                const labelListList& addr,
                scalarListList& wght,
                scalarList& symWght,
                scalarField& wghtSum,
                const bool conformal,
                const bool output,
                const scalar lowWeightTol,
                const scalar maxScale
            );

            void cropWeights
            (
                const scalarField& patchAreas,
                const word& patchName,
                scalarListList& wght,
                scalarField& wghtSum
            );


        // Constructor helpers

            static void agglomerate
            (
                const autoPtr<mapDistribute>& targetMap,
                const scalarField& fineSrcMagSf,
                const labelListList& fineSrcAddress,
                const scalarListList& fineSrcWeights,
                const scalarList& fineSrcSymWeights,

                const labelList& sourceRestrictAddressing,
                const labelList& targetRestrictAddressing,
                const scalar maxScale,

                scalarField& srcMagSf,
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarList& srcSymWeights,
                scalarField& srcWeightsSum,
                autoPtr<mapDistribute>& tgtMap
            );

            void constructFromSurface
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr
            );

            //- Construct overlap and nonOverlap addressing lists
            void constructOverAndNonOverlappingAddressing();

            //- Construct overlap and nonOverlap addressing lists
            //  Checks the sizes of provided src and tgtPatch
            void constructOverAndNonOverlappingAddressing
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            );

public:

    // Constructors

        //- Construct from components
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const scalar& relTol_ = relTolDef_,
            const scalar& tol_ = tolDef_,
            const scalar& cosMatchAngle = cosMatchAngleDef_,
            const scalar& maxAMIWeightScale = maxAMIWeightScaleDef_,
            const bool output = outputDef_
        );

            void visualiseWeights
            (
                const SourcePatch&,
                const TargetPatch&,
                const bool,
                const word&
            ) const;

        //- Construct from components
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const word& methodName =
                interpolationMethodToWord(imFaceAreaWeight),
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const scalar& relTol_ = relTolDef_,
            const scalar& tol_ = tolDef_,
            const scalar& cosMatchAngle = cosMatchAngleDef_,
            const scalar& maxAMIWeightScale = maxAMIWeightScaleDef_,
            const bool output = outputDef_
        );

        //- Construct from components, with projection surface
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const autoPtr<searchableSurface>& surf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const scalar& relTol_ = relTolDef_,
            const scalar& tol_ = tolDef_,
            const scalar& cosMatchAngle = cosMatchAngleDef_,
            const scalar& maxAMIWeightScale = maxAMIWeightScaleDef_,
            const bool output = outputDef_
        );

        //- Construct from components, with projection surface
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const autoPtr<searchableSurface>& surf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const word& methodName =
                interpolationMethodToWord(imFaceAreaWeight),
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const scalar& relTol_ = relTolDef_,
            const scalar& tol_ = tolDef_,
            const scalar& cosMatchAngle = cosMatchAngleDef_,
            const scalar& maxAMIWeightScale = maxAMIWeightScaleDef_,
            const bool output = outputDef_
        );

        //- Construct from agglomeration of AMIInterpolation. Agglomeration
        //  passed in as new coarse size and addressing from fine from coarse
        AMIInterpolation
        (
            const AMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
            const labelList& sourceRestrictAddressing,
            const labelList& neighbourRestrictAddressing
        );


    //- Destructor
    ~AMIInterpolation();

    // Typedef to SourcePatch type this AMIInterplation is instantiated on
    typedef SourcePatch sourcePatchType;

    // Typedef to TargetPatch type this AMIInterplation is instantiated on
    typedef TargetPatch targetPatchType;


    // Member Functions

        // Access

            //- Access to the up-to-date flag
            inline bool upToDate() const;

            //- Access to the up-to-date flag
            inline bool& upToDate();

            //- Set to -1, or the processor holding all faces (both sides) of
            //  the AMI
            inline label singlePatchProc() const;

            //- Threshold weight below which interpolation is deactivated
            inline scalar lowWeightCorrection() const;

            //- Return true if employing a 'lowWeightCorrection'
            inline bool applyLowWeightCorrection() const;


            // Source patch

                //- Return const access to source patch face areas
                inline const scalarField& srcMagSf() const;

                //- Return access to source patch face areas
                inline scalarField& srcMagSf();

                //- Return const access to source patch face normals
                inline const vectorField& srcNf() const;

                //- Return const access to source patch addressing
                inline const labelListList& srcAddress() const;

                //- Return access to source patch addressing
                inline labelListList& srcAddress();

                //- Return const access to source patch weights
                inline const scalarListList& srcWeights() const;

                //- Return access to source patch weights
                inline scalarListList& srcWeights();

                //- Return const access to source patch weights
                inline const scalarList& srcSymWeights() const;

                //- Return access to source patch weights
                inline scalarList& srcSymWeights();

                //- Return const access to normalisation factor of source
                //  patch weights (i.e. the sum before normalisation)
                inline const scalarField& srcWeightsSum() const;

                //- Return access to normalisation factor of source
                //  patch weights (i.e. the sum before normalisation)
                inline scalarField& srcWeightsSum();

                //- Return const access to source patch face centroids
                inline const pointListList& srcCentroids() const;

                //- Return access to source patch face centroids
                inline pointListList& srcCentroids();

                //- Source map pointer - valid only if singlePatchProc = -1
                //  This gets source data into a form to be consumed by
                //  tgtAddress, tgtWeights
                inline const mapDistribute& srcMap() const;


            // Target patch

                //- Return const access to target patch face areas
                inline const scalarField& tgtMagSf() const;

                //- Return access to target patch face areas
                inline scalarField& tgtMagSf();

                //- Return const access to target patch face normals
                inline const vectorField& tgtNf() const;

                //- Return const access to target patch addressing
                inline const labelListList& tgtAddress() const;

                //- Return access to target patch addressing
                inline labelListList& tgtAddress();

                //- Return const access to target patch weights
                inline const scalarListList& tgtWeights() const;

                //- Return access to target patch weights
                inline scalarListList& tgtWeights();

                //- Return const access to source patch weights
                inline const scalarList& tgtSymWeights() const;

                //- Return access to source patch weights
                inline scalarList& tgtSymWeights();

                //- Return const access to normalisation factor of target
                //  patch weights (i.e. the sum before normalisation)
                inline const scalarField& tgtWeightsSum() const;

                //- Return access to normalisation factor of target
                //  patch weights (i.e. the sum before normalisation)
                inline scalarField& tgtWeightsSum();

                //- Target map pointer -  valid only if singlePatchProc=-1.
                //  This gets target data into a form to be consumed by
                //  srcAddress, srcWeights
                inline const mapDistribute& tgtMap() const;


                //- Coupled solver inline access functions

                 inline const labelList& overlapTargetFaces() const;

                 inline const labelList& nonOverlapTargetFaces() const;

                 inline const labelList& overlapSourceFaces() const;

                 inline const labelList& nonOverlapSourceFaces() const;


        // Manipulation

            //- Update addressing and weights
            void update
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch
            );


            //- Set the maps, addresses and weights from an external source
            void reset
            (
                autoPtr<mapDistribute>&& srcToTgtMap,
                autoPtr<mapDistribute>&& tgtToSrcMap,
                labelListList&& srcAddress,
                scalarListList&& srcWeights,
                labelListList&& tgtAddress,
                scalarListList&& tgtWeights
            );


            //- Append additional addressing and weights
            void append
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const label tranfI = 0
            );

            //- Normalise the weights
            void normaliseWeights(const bool conformal, const bool output);

            //- Crom weights more than one
            void cropWeights();

            void setTransformations
            (
                const List<transformer>&,
                const List<transformer>&
            );


        // Evaluation

            // Low-level

                //- Interpolate from target to source with supplied op
                //  to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToSource
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;

                //- Interpolate from source to target with supplied op
                //  to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToTarget
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;


            //- Interpolate from target to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target to source
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

        // Point intersections

            //- Return source patch face index of point on target patch face
            label srcPointFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const vector& n,
                const label tgtFacei,
                point& tgtPoint
            )
            const;

            //- Return target patch face index of point on source patch face
            label tgtPointFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const vector& n,
                const label srcFacei,
                point& srcPoint
            )
            const;


        // Checks

            //- Write face connectivity as OBJ file
            void writeFaceConnectivity
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const labelListList& srcAddress
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolation/AMIInterpolation/AMIInterpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIInterpolation/AMIInterpolation/AMIInterpolation.C"
    #include "AMIInterpolation/AMIInterpolation/AMIInterpolationParallelOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
