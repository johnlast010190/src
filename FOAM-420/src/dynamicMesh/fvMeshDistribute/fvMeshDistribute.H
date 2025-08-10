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
    (c) 2011-2018 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2022-2023 Esi Ltd.

Class
    Foam::fvMeshDistribute

Description
    Send/receive parts of mesh + fvFields to neighbouring processors.
    Used in load balancing and mesh redistribution.

    Input is per local cell the processor it should move to. Moves meshes
    and volFields/surfaceFields, returning a map which can be used to
    distribute other.

    Notes:
    - does not handle cyclics. Will probably handle separated proc patches.
    - if all cells move off processor also all its processor patches will
      get deleted so comms might be screwed up (since e.g. globalMeshData
      expects procPatches on all).
    - initial mesh has to have procPatches last and all normal patches common
      to all processors and in the same order. This is checked.
    - faces are matched topologically but points on the faces are not. So
      expect problems -on separated patches (cyclics?) -on zero sized processor
      edges.
    - handles non-conformal cyclic patches.

SourceFiles
    fvMeshDistribute.C
    fvMeshDistributeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistribute_H
#define fvMeshDistribute_H

#include "fields/Fields/Field/Field.H"
#include "fvMeshSubset/fvMeshSubset.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "VectorN/finiteVolume/fields/surfaceFields/surfaceVectorNFields.H"
#include "VectorN/DimensonedFields/DimensionedDiagTensorNFields.H"
#include "VectorN/DimensonedFields/DimensionedSphericalTensorNFields.H"
#include "VectorN/DimensonedFields/DimensionedTensorNFields.H"
#include "VectorN/DimensonedFields/DimensionedVectorNFields.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapAddedPolyMesh;

// default parameters for distribute
static PtrList<fvMesh> dummyMeshList;
static PtrList<mapDistributePolyMesh> dummyMapList;

/*---------------------------------------------------------------------------*\
                      Class fvMeshDistribute Declaration
\*---------------------------------------------------------------------------*/

class fvMeshDistribute
{
    // Private data

        //- Underlying fvMesh
        fvMesh& mesh_;

        //- Absolute merging tolerance (constructing meshes gets done using
        //  geometric matching)
        const scalar mergeTol_;

        //- Total number of partitions after distribution
        label nFinalProcs_;


    // Private Member Functions

        static void inplaceRenumberWithFlip
        (
            const labelUList& oldToNew,
            const bool oldToNewHasFlip,
            const bool lstHasFlip,
            labelUList& lst
        );

        //- Find indices with value
        static labelList select
        (
            const bool selectEqual,
            const labelList& values,
            const label value
        );

        //- Check all procs have same names and in exactly same order.
        static void checkEqualWordList(const string&, const wordList&);


        // Patch handling

            //- Find patch to put exposed faces into.
            label findNonEmptyPatch() const;

            //- Save boundary fields
            template<class T, class Mesh>
            void saveBoundaryFields
            (
                PtrList<FieldField<fvsPatchField, T>>& bflds,
                const fvMesh& mesh
            ) const;

            //- Map boundary fields
            template<class T, class Mesh>
            void mapBoundaryFields
            (
                const mapPolyMesh& map,
                const PtrList<FieldField<fvsPatchField, T>>& oldBflds,
                fvMesh& mesh
            );

            //- Save internal fields of surfaceFields
            template<class T>
            void saveInternalFields(PtrList<Field<T>>& iflds) const;

            //- Set value of patch faces resulting from internal faces
            template<class T>
            void mapExposedFaces
            (
                const mapPolyMesh& map,
                const PtrList<Field<T>>& oldFlds
            );

            //- Init patch fields of certain type
            template<class GeoField, class PatchFieldType>
            void initPatchFields
            (
                const typename GeoField::value_type& initVal,
                fvMesh& mesh
            );

            //- Call correctBoundaryConditions on fields
            template<class GeoField>
            void correctBoundaryConditions();

            //- Reset fvPatchfield updated_ flag for fields
            template<class GeoField>
            void resetUpdate();

            //- Delete all processor patches. Move any processor faces into
            //  patchi.
            autoPtr<mapPolyMesh> deleteProcPatches(const label patchi);

            //- Repatch the mesh. This is only necessary for the proc
            //  boundary faces. newPatchID is over all boundary faces: -1 or
            //  new patchID. constructFaceMap is being adapted for the
            //  possible new face position (since proc faces get automatically
            //  matched)
            autoPtr<mapPolyMesh> repatch
            (
                const labelList& newPatchID,
                fvMesh& mesh,
                labelListList& constructFaceMap
            );

            //- Merge any local points that were remotely coupled.
            //  constructPointMap is adapted for the new point labels.
            autoPtr<mapPolyMesh> mergeSharedPoints
            (
                const labelList& pointToGlobalMaster,
                fvMesh& mesh,
                labelListList& constructPointMap
            );


        // Coupling information

            //- Construct the local environment of all boundary faces.
            void getCouplingData
            (
                const labelList& distribution,
                labelList& sourceFace,
                labelList& sourceProc,
                labelList& sourcePatch,
                labelList& sourceNbrPatch,
                labelList& sourceNewProc,
                labelList& sourcePointMaster
            ) const;

            // Subset the neighbourCell/neighbourProc fields
            static void subsetCouplingData
            (
                const fvMesh& mesh,
                const labelList& pointMap,
                const labelList& faceMap,
                const labelList& cellMap,

                const labelList& oldDistribution,
                const labelList& oldFaceOwner,
                const labelList& oldFaceNeighbour,
                const label oldInternalFaces,

                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,
                const labelList& sourceNbrPatch,
                const labelList& sourceNewProc,
                const labelList& sourcePointMaster,

                labelList& subFace,
                labelList& subProc,
                labelList& subPatch,
                labelList& subNbrPatch,
                labelList& subNewProc,
                labelList& subPointMaster
            );

            //- Find cells on mesh whose faceID/procID match the neighbour
            //  cell/proc of domainMesh. Store the matching face.
            static void findCouples
            (
                const primitiveMesh&,
                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,

                const label domain,
                const primitiveMesh& domainMesh,
                const labelList& domainFace,
                const labelList& domainProc,
                const labelList& domainPatch,

                labelList& masterCoupledFaces,
                labelList& slaveCoupledFaces
            );

            //- Map data on boundary faces to new mesh (resulting from adding
            //  two meshes)
            static labelList mapBoundaryData
            (
                const primitiveMesh& mesh,      // mesh after adding
                const mapAddedPolyMesh& map,
                const labelList& boundaryData0, // mesh before adding
                const label nInternalFaces1,
                const labelList& boundaryData1  // added mesh
            );

            //- Map point data to new mesh (resulting from adding two meshes)
            static labelList mapPointData
            (
                const primitiveMesh& mesh,      // mesh after adding
                const mapAddedPolyMesh& map,
                const labelList& boundaryData0, // on mesh before adding
                const labelList& boundaryData1  // on added mesh
            );

        // Other

            //- Remove cells. Add all exposed faces to patch oldInternalPatchi
            autoPtr<mapPolyMesh> doRemoveCells
            (
                const labelList& cellsToRemove,
                const label oldInternalPatchi
            );

            //- Add processor patches. Changes mesh and returns per neighbour
            //  proc the processor patchID.
            void addProcPatches
            (
                const labelList& nbrProc,         // Processor that neighbour
                                                  // is now on
                const labelList& referPatchID,    // Original patch ID (or -1)
                const labelList& referNbrPatchID, // Original neighbour patch
                                                  // ID (or -1)
                const label myMeshNo,
                fvMesh& mesh,
                List<Map<label>>& procPatchID
            );

            //- Add non-conformal processor cyclic patches. Changes the mesh.
            void addNccProcPatches
            (
                PtrList<fvMesh>& myProcMeshes,
                const List<labelList>& meshNoInProc
            );

            //- Construct poly-face addressing between original and final
            //  processor meshes, for a given non-conformal cyclic patch.
            List<List<DynamicList<label>>> procPolyFaceAddr
            (
                const List<labelList>& subMap,
                const surfaceLabelField::Boundary origPolyFacesBf,
                const labelList& origFaceOwner,
                const labelList& cellToProc,
                const label nccPatchID,
                const label nccNbrPatchID
            ) const;

            //- If the original processor meshes are non-conformal, then map
            //  the non-conformal addressing across and unconform all the
            //  redistributed processor meshes.
            void unconform
            (
                PtrList<fvMesh>& myProcMeshes,
                const List<labelList>& meshNoInProc,
                const PtrList<mapDistributePolyMesh>& distMaps,
                const surfaceLabelField::Boundary origPolyFacesBf,
                const labelList& origFaceOwner,
                const labelList& cellToProc
            );

            //- Get boundary faces to be repatched. Is -1 or new patchID
            static labelList getBoundaryPatch
            (
                const labelList& neighbourNewProc,  // new processor per b. face
                const labelList& referPatchID,      // -1 or original patch
                const List<Map<label>>& procPatchID,// patchID
                const label myMeshNo
            );

            //- Send mesh and coupling data.
            static void sendMesh
            (
                const label domain,
                const fvMesh& mesh,
                const wordList& pointZoneNames,
                const wordList& facesZoneNames,
                const wordList& cellZoneNames,
                const labelList& sourceFace,
                const labelList& sourceProc,
                const labelList& sourcePatch,
                const labelList& sourceNbrPatch,
                const labelList& sourceNewProc,
                const labelList& sourcePointMaster,
                Ostream& toDomain
            );
            //- Send subset of fields
            template<class GeoField>
            static void sendFields
            (
                const label domain,
                const wordList& fieldNames,
                const fvMeshSubset&,
                Ostream& toNbr
            );

            //- Receive mesh. Opposite of sendMesh
            static autoPtr<fvMesh> receiveMesh
            (
                const label domain,
                const wordList& pointZoneNames,
                const wordList& facesZoneNames,
                const wordList& cellZoneNames,
                const Time& runTime,
                labelList& domainSourceFace,
                labelList& domainSourceProc,
                labelList& domainSourcePatch,
                labelList& domainSourceNbrPatch,
                labelList& domainSourceNewProc,
                labelList& domainSourcePointMaster,
                Istream& fromNbr
            );

            //- Initialize empty fields
            template<class GeoField>
            void initFields
            (
                fvMesh& baseMesh,
                fvMesh& mesh,
                const autoPtr<fvMeshSubset>& subsetterPtr,
                PtrList<GeoField>& fields
            );

            //- Receive fields. Opposite of sendFields
            template<class GeoField>
            static void receiveFields
            (
                const label domain,
                const wordList& fieldNames,
                fvMesh&,
                PtrList<GeoField>&,
                const dictionary& fieldDicts
            );

            //- Helper function: get the number of vertices face f needs to be
            //  rotated such that its f[0] point aligns with given anchor
            //  (within a specified tolerance)
            label getRotation
            (
                const pointField& points,
                const face& f,
                const point& anchor,
                const scalar tol
            );

            //- Reorder processor faces
            void reorderCoupledFaces
            (
                PtrList<fvMesh>& myProcMeshes,
                const labelList& meshToProc,
                const labelListList& meshNoInProc
            );

            template<typename... Args>
            void receiveAllFields
            (
                const label sendProc,
                fvMesh& domainMesh,
                const dictionary& fieldDicts,

                Args&... args
            );

            template<typename... Args, size_t... Is>
            void receiveAllFieldsImpl
            (
                const label sendProc,
                fvMesh& domainMesh,
                const dictionary& fieldDicts,
                std::index_sequence<Is...>,

                Args&... args
            );

            //- No copy construct
            fvMeshDistribute(const fvMeshDistribute&) = delete;

            //- No copy assignment
            void operator=(const fvMeshDistribute&) = delete;

public:

    ClassName("fvMeshDistribute");


    // Constructors

        //- Construct from mesh and absolute merge tolerance
        fvMeshDistribute(fvMesh& mesh, const scalar mergeTol);


    // Member Functions

        //- Helper function: count cells per processor in wanted distribution
        static labelList countCells(const labelList&);
        static labelList countCells(const labelList&, const label);

        //- Merge wordlists over all processors
        static wordList mergeWordList(const wordList&);

        //- Calculates which mesh goes to which processor
        static labelListList calcMeshToProcMap
        (
            const label nFinalProcs,
            labelList& nMeshesInProc,
            labelList& meshToProc
        );

        //- Set the final number of partitions
        void setNFinalProcs(const label nDestProcs);

        //- Initialize list of empty meshes and fields for each MPI rank
        // when distributing to higher number of processors
        void initMeshList
        (
            const fileName& meshSubDir,
            const label nDestProcs,
            fvMesh& baseMesh,

            PtrList<Time>&                 myProcTimes,
            PtrList<fvMesh>&               myProcMeshes,
            List<PtrList<volScalarField>>& volScalarFields,
            List<PtrList<volVectorField>>& volVectorFields
        );

        //- Send cells to neighbours according to distribution
        //  (for every cell the new proc)
        autoPtr<mapDistributePolyMesh> distribute
        (
            const labelList& dist,
            PtrList<fvMesh>& myProcMeshes = dummyMeshList,
            PtrList<mapDistributePolyMesh>& mapsDist = dummyMapList
        );

        // Debugging

            //- Print some info on coupling data
            static void printCoupleInfo
            (
                const primitiveMesh&,
                const labelList&,
                const labelList&,
                const labelList&,
                const labelList&
            );

            //- Print some field info
            template<class GeoField>
            static void printFieldInfo(const fvMesh&);

            //- Print some info on mesh.
            static void printMeshInfo(const fvMesh&);

            //- Generate a test field on faces
            static tmp<surfaceScalarField> generateTestField(const fvMesh&);

            //- Check whether field consistent with face orientation
            static void testField(const surfaceScalarField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshDistribute/fvMeshDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
