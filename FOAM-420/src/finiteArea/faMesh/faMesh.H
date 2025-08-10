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
    (c) 2016-2022 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

Class
    Foam::faMesh

Description
    Finite area mesh.  Used for 2-D non-Euclidian finite area method.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faMesh_H
#define faMesh_H

#include "meshes/GeoMesh/GeoMesh.H"
#include "meshes/MeshObject/MeshObject.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/lduMesh/lduMesh.H"
#include "faMesh/faBoundaryMesh/faBoundaryMesh.H"
#include "meshes/meshShapes/edge/edgeList.H"
#include "meshes/meshShapes/face/faceList.H"
#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "interpolation/edgeInterpolation/edgeInterpolation.H"
#include "primitives/ints/lists/labelIOList.H"
#include "fields/Fields/scalarField/scalarIOField.H"
#include "fields/FieldFields/FieldField/FieldFields.H"
#include "faMesh/faGlobalMeshData/faGlobalMeshData.H"
#include "finiteArea/faSchemes/faSchemes.H"
#include "faSolution/faSolution.H"
#include "meshes/data/data.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faMeshLduAddressing;
class faMeshMapper;

/*---------------------------------------------------------------------------*\
                             Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh
:
    public GeoMesh<polyMesh>,
    public MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>,
    public lduMesh,
    public edgeInterpolation,
    public data
{
    // Private data

        //- Face labels
        labelIOList faceLabels_;

        //- Boundary mesh
        faBoundaryMesh boundary_;


        // Primitive mesh data

            //- Edges, addressing into local point list
            edgeList edges_;

            //- Edge owner
            labelList edgeOwner_;

            //- Edge neighbour
            labelList edgeNeighbour_;


        // Primitive size data

            //- Number of points
            mutable label nPoints_;

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal edges
            mutable label nInternalEdges_;

            //- Number of faces
            mutable label nFaces_;


        // Communication support

            //- Communicator used for parallel communication
            label comm_;


        // When constructing from a single patch, stores a reference to a
        // sub-registry to allow multiple faMeshes on a single polyMesh
        const objectRegistry& registry_;

    // Demand-driven data

        //- Primitive patch
        mutable indirectPrimitivePatch* patchPtr_;

        //- Ldu addressing data
        mutable faMeshLduAddressing* lduPtr_;

        //- Current time index for motion
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Face areas
        mutable DimensionedField<scalar, areaMesh>* SPtr_;

        //- Face areas old time level
        mutable DimensionedField<scalar, areaMesh>* S0Ptr_;

        //- Face areas old-old time level
        mutable DimensionedField<scalar, areaMesh>* S00Ptr_;

        //- Patch starts in the edge list
        mutable labelList* patchStartsPtr_;

        //- Edge length vectors
        mutable edgeVectorField* LePtr_;

        //- Mag edge length vectors
        mutable edgeScalarField* magLePtr_;

        //- Face centres
        mutable areaVectorField* centresPtr_;

        //- Edge centres
        mutable edgeVectorField* edgeCentresPtr_;

        //- Face area normals
        mutable areaVectorField* faceAreaNormalsPtr_;

        //- Edge area normals
        mutable edgeVectorField* edgeAreaNormalsPtr_;

        //- Edge area normals
        mutable vectorField* pointAreaNormalsPtr_;

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;

        //- Edge transformation tensors
        mutable FieldField<Field, tensor>* edgeTransformTensorsPtr_;

        //- Whether point normals must be corrected for a patch
        mutable boolList* correctPatchPointNormalsPtr_;

        //- Face motion fluxes
        mutable edgeScalarField* phiPtr_;

        //- faSchemes created on demand
        mutable autoPtr<faSchemes> faSchemes_;

        //- faSolution created on demand
        mutable autoPtr<faSolution> faSolution_;


        // Other mesh-related data

            //- Parallel info
            mutable faGlobalMeshData* globalMeshDataPtr_;

    // Static Private Data

        //- Use quadrics fit
        static const int quadricsFit_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faMesh(const faMesh&);

        //- Disallow default bitwise assignment
        void operator=(const faMesh&);


        //- Set primitive mesh data
        void setPrimitiveMeshData();

        // Private member functions to calculate demand driven data

            //- Calculate ldu addressing
            void calcLduAddressing() const;

            //- Calculate patch starts in the edge list
            void calcPatchStarts() const;

            //- Calculate edge lengths
            void calcLe() const;

            //- Calculate mag edge lengths
            void calcMagLe() const;

            //- Calculate face centres
            void calcAreaCentres() const;

            //- Calculate edge centres
            void calcEdgeCentres() const;

            //- Calculate face areas
            void calcS() const;

            //- Calculate face area normals
            void calcFaceAreaNormals() const;

            //- Calculate edge area normals
            void calcEdgeAreaNormals() const;

            //- Calculate point area normals
            void calcPointAreaNormals() const;

            //- Calculate point area normals by quadrics fit
            void calcPointAreaNormalsByQuadricsFit() const;

            //- Calculate face curvatures
            void calcFaceCurvatures() const;

            //- Calculate edge transformation tensors
            void calcEdgeTransformTensors() const;

            //- Clear geometry but not the face areas
            void clearGeomNotOldAreas() const;

            //- Clear geometry like clearGeomNotOldAreas but recreate any
            //  geometric demand-driven data that was set
            void updateGeomNotOldAreas();

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

            //- Preserve old volume(s)
            void storeOldSurf(const scalarField&);

            //- Clear demand-driven data
            void clearOut() const;

            tmp<scalarField> sweptAreas
            (
                const pointField& newPoints,
                const pointField& oldPoints
            );

public:

    // Public typedefs

        typedef faMesh Mesh;
        typedef faBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("faMesh");


    //- Return the mesh sub-directory name (usually "faMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from polyMesh
        explicit faMesh(const polyMesh& m);

        //- Construct from components without boundary.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const polyMesh& m,
            const labelList& faceLabels
        );

        //- Construct from finite area mesh definition file
        faMesh
        (
            const polyMesh& m,
            const fileName& defFile
        );

        //- Construct from polyPatch
        faMesh
        (
            const polyMesh& m,
            const label polyPatchID
        );

        faMesh
        (
            const polyMesh& pMesh,
            const word& patchName
        );

        faMesh
        (
            const polyMesh& pMesh,
            const wordReList& patchNames
        );


    // Destructor

        virtual ~faMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addFaPatches(const List<faPatch*> &);


        // Database

            //- Return access to polyMesh
            const polyMesh& mesh() const
            {
                return MeshObject<polyMesh, UpdateableMeshObject, faMesh>::mesh();
            }

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return reference to time
            const Time& time() const;

            //- Return the current instance directory for points
            //  Used in the consruction of gemometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;


            //- Mesh size parameters

                inline label nPoints() const
                {
                    return nPoints_;
                }

                inline label nEdges() const
                {
                    return nEdges_;
                }

                inline label nInternalEdges() const
                {
                    return nInternalEdges_;
                }

                inline label nFaces() const
                {
                    return nFaces_;
                }

            // Primitive mesh data

                //- Return mesh points
                const pointField& points() const;

                //- Return edges
                const edgeList& edges() const;

                //- Return faces
                const faceList& faces() const;

                //- Edge owner addresing
                inline const labelList& edgeOwner() const
                {
                    return edgeOwner_;
                }

                //- Edge neighbour addressing
                inline const labelList& edgeNeighbour() const
                {
                    return edgeNeighbour_;
                }


        // Access

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const override
            {
                return true;
            }

            //- Return reference to the mesh database
            const objectRegistry& thisDb() const override;

            //- Return reference to name
            //  Note: name() is currently ambiguous due to derivation from
            //  surfaceInterpolation
            const word& name() const
            {
                //Dont know how to access the name
                return mesh().name();
            }

            //- Return constant reference to boundary mesh
            const faBoundaryMesh& boundary() const;

            //- Return faMesh face labels
            const labelList& faceLabels() const
            {
                return faceLabels_;
            }


            //- Return parallel info
            const faGlobalMeshData& globalData() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const override;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const override
            {
                return boundary().interfaces();
            }

//            //- Return communicator used for parallel communication
            virtual label comm() const override
            {
                return mesh().comm();
            }


            //- Internal face owner
            const labelUList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return true if given edge label is internal to the mesh
            inline bool isInternalEdge(const label edgeIndex) const
            {
                return edgeIndex < nInternalEdges();
            }


        // Mesh motion and mophing

            //- Is mesh moving
            bool moving() const
            {
                return mesh().moving();
            }

            //- Update after mesh motion
            bool movePoints() override;

            //- Update after topo change
            void updateMesh(const mapPolyMesh&) override;


        // Mapping

            //- Map all fields in time using given map.
            virtual void mapFields(const faMeshMapper& mapper) const;

            //- Map face areas in time using given map.
            virtual void mapOldAreas(const faMeshMapper& mapper) const;


        // Demand-driven data

            //- Return constant reference to primitive patch
            const indirectPrimitivePatch& patch() const;

            //- Return reference to primitive patch
            indirectPrimitivePatch& patch();

            //- Return patch starts
            const labelList& patchStarts() const;

            //- Return edge length vectors
            const edgeVectorField& Le() const;

            //- Return edge length magnitudes
            const edgeScalarField& magLe() const;

            //- Return face centres as areaVectorField
            const areaVectorField& areaCentres() const;

            //- Return edge centres as edgeVectorField
            const edgeVectorField& edgeCentres() const;

            //- Return face areas
            virtual
            const DimensionedField<scalar, areaMesh>& S() const;

            //- Return old-time face areas
            virtual
            const DimensionedField<scalar, areaMesh>& S0() const;

            //- Return old-old-time face areas
            virtual
            const DimensionedField<scalar, areaMesh>& S00() const;

            //- Return sub-cycle face areas
            virtual
            tmp<DimensionedField<scalar, areaMesh>> Ssc() const;

            //- Return sub-cycl old-time face areas
            virtual
            tmp<DimensionedField<scalar, areaMesh>> Ssc0() const;

            //- Return the vector of solved-for directions in mesh.
            //- Calculating through polyMesh class
            const Vector<label>& solutionD() const;

            //- Return face area normals
            const areaVectorField& faceAreaNormals() const;

            //- Return edge area normals
            const edgeVectorField& edgeAreaNormals() const;

            //- Return point area normals
            const vectorField& pointAreaNormals() const;

            //- Return face curvatures
            const areaScalarField& faceCurvatures() const;

            //- Return edge transformation tensors
            const FieldField<Field, tensor>& edgeTransformTensors() const;

            //- Return internal point labels
            labelList internalPoints() const;

            //- Return boundary point labels
            labelList boundaryPoints() const;

            //- Return edge length correction
            tmp<edgeScalarField> edgeLengthCorrection() const;

            //- Whether point normals should be corrected for a patch
            bool correctPatchPointNormals(const label patchID) const;

            //- Set whether point normals should be corrected for a patch
            boolList& correctPatchPointNormals() const;

            //- Return face edge motion fluxes
            const edgeScalarField& phi() const;

            //- Clear mesh motion flux
            void clearMeshPhi();

            //- Return the fvSchemes
            const faSchemes& schemes() const;

            //- Return the fvSchemes
            const faSolution& solution() const;

            //- Return the fvSchemes
            faSchemes& schemes();

            //- Return the fvSchemes
            faSolution& solution();

        //- Write mesh
        virtual bool write(const bool valid = true) const override;

        // Member Operators

//            //- Return reference to polyMesh
//            const polyMesh& operator()() const
//            {
//                return MeshObject<polyMesh, UpdateableMeshObject, faMesh>::mesh();
//            }

    // Member Operators

        bool operator!=(const faMesh& m) const;

        bool operator==(const faMesh& m) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faMesh/faPatches/faPatch/faPatchFaMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
