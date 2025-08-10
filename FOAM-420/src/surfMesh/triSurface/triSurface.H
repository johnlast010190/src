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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

Class
    Foam::triSurface

Description
    Triangulated surface description with patch information.

SourceFiles
    triSurface.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface_H
#define triSurface_H

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "primitives/bools/lists/boolList.H"
#include "triSurface/patches/geometricSurfacePatchList.H"
#include "triSurface/patches/surfacePatchList.H"
#include "meshes/meshShapes/triFace/triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class IFstream;
class surfZone;

// Forward declaration of friend functions and operators

class triSurface;
Istream& operator>>(Istream&, triSurface&);
Ostream& operator<<(Ostream&, const triSurface&);


/*---------------------------------------------------------------------------*\
                           Class triSurface Declaration
\*---------------------------------------------------------------------------*/

class triSurface
:
    public PrimitivePatch<::Foam::List<labelledTri>, pointField>
{
    // Private typedefs

    //- Typedefs for convenience
        typedef labelledTri Face;
        typedef PrimitivePatch<::Foam::List<labelledTri>, pointField>
            ParentType;


    // Private data

        //- The number of bytes in the STL header
        static const int STLheaderSize = 80;

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        geometricSurfacePatchList patches_;

        static const wordHashSet readTypes_;
        static const wordHashSet writeTypes_;


    // Demand driven private data.

        //- Edge-face addressing (sorted)
        mutable labelListList* sortedEdgeFacesPtr_;

        //- Label of face that 'owns' edge (i.e. e.vec() is righthanded walk
        //  along face)
        mutable labelList* edgeOwnerPtr_;

        mutable pointField* points0Ptr_;


    // Private Member Functions

        //- Calculate sorted edgeFaces
        void calcSortedEdgeFaces() const;

        //- Calculate owner
        void calcEdgeOwner() const;

        //- Sets default values for patches
        void setDefaultPatches();

        //- Apply dictionary transforms
        void doTransforms
        (
            const word& name,
            const dictionary& dict
        );

        //- Function to stitch the triangles by removing duplicate points.
        //  Returns true if any points merged
        bool stitchTriangles
        (
            const scalar tol = SMALL,
            const bool verbose = false
        );

        //- Function to remove unused points.
        bool removeUnusedPoints
        (
            const bool verbose = false
        );

        //- Read in Foam format
        bool read(Istream&);

        //- Generic read routine. Chooses reader based on extension.
        bool read
        (
            const fileName&,
            const word& ext,
            const dictionary& dict,
            const bool check = true
        );

        bool readSTL(const fileName&, bool forceBinary=false);
        bool readGTS(const fileName&);
        bool readOBJ(const fileName&);
        bool readOFF(const fileName&);
        bool readTRI(const fileName&);
        bool readAC(const fileName&);
        bool readNAS(const fileName&);
        bool readINP(const fileName&);
        bool readNTL(const fileName&);
        bool readVTK(const fileName&);
        bool readBTS(const fileName&);
        bool readALEMBIC(const fileName&, const dictionary& dict);

        //- Generic write routine. Chooses writer based on extension.
        void write
        (
            const fileName&,
            const word& ext,
            const bool sort,
            const bool keepEmpty
        ) const;

        //- Write to Ostream in ASCII STL format.
        //  Each region becomes 'solid' 'endsolid' block.
        void writeSTLASCII
        (
               const bool writeSorted,
               const bool keepEmpty,
               Ostream&
        ) const;

        //- Write to std::ostream in BINARY STL format
        void writeSTLBINARY(std::ostream&) const;

        //- Write to Ostream in ASCII NAS format.
        void writeNASASCII(Ostream&) const;

        //- Write to Ostream in GTS (Gnu Tri Surface library)
        //  format.
        void writeGTS(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in OBJ (Lightwave) format.
        //  writeSorted=true: sort faces acc. to region and write as single
        //  group. =false: write in normal order.
        void writeOBJ(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in OFF (Geomview) format.
        //  writeSorted=true: sort faces acc. to region and write as single
        //  group. =false: write in normal order.
        void writeOFF(const bool writeSorted, Ostream&) const;

        //- Write to VTK legacy format.
        void writeVTK(const bool writeSorted, std::ostream& os) const;

        //- Write to Ostream in TRI (AC3D) format
        //  Ac3d .tri format (unmerged triangle format)
        void writeTRI(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in SMESH (tetgen) format
        void writeSMESH(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in AC3D format. Always sorted by patch.
        void writeAC(Ostream&) const;

        //- Write to Ostream in ASCII NAS format.
        void writeINP(Ostream&) const;


    // Static private functions

        //- Convert faces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const faceList&,
            const label defaultRegion = 0
        );

        //- Convert triFaces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const triFaceList&,
            const label defaultRegion = 0
        );

        //- Helper function to print triangle info
        static void printTriangle
        (
            Ostream&,
            const Foam::string& pre,
            const labelledTri&,
            const pointField&
        );

        //- Read non-comment line
        static string getLineNoComment(IFstream&);


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face>&>(*this);
        }


public:

    // Public typedefs

        //- Placeholder only, but do not remove - it is needed for GeoMesh
        typedef bool BoundaryMesh;


        //- Runtime type information
        ClassName("triSurface");


    // Static

        //- Name of triSurface directory to use.
        static fileName triSurfInstance(const Time&);

        //- Can we read this file format?
        static bool canRead(const fileName& name, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static const wordHashSet& readTypes();
        static const wordHashSet& writeTypes();


    // Constructors

        //- Construct null
        triSurface();

        //- Construct from triangles, patches, points.
        triSurface
        (
            const List<labelledTri>&,
            const geometricSurfacePatchList&,
            const pointField&
        );

        //- Construct from triangles, patches, points. Reuse storage.
        triSurface
        (
            List<labelledTri>&,
            const geometricSurfacePatchList&,
            pointField&,
            const bool reuse
        );

        //- Construct by transferring (triangles, points) components.
        triSurface
        (
            const Xfer<List<labelledTri>>&,
            const geometricSurfacePatchList&,
            const Xfer<List<point>>&
        );

        //- Construct from triangles, points. Set patch names to default.
        triSurface(const List<labelledTri>&, const pointField&);

        //- Construct from triangles, points. Set region to 0 and default
        //  patchName.
        triSurface(const triFaceList&, const pointField&);

        //- Construct from file name (uses extension to determine type)
        //  Optional (positive, non-zero) point scaling is possible.
        triSurface
        (
            const fileName&,
            const scalar scaleFactor = -1,
            const dictionary dict = dictionary()
        );

        //- Construct and apply local and global transforms
        triSurface
        (
            const IOobject& io,
            const fileName& fname,
            dictionary dict = dictionary()
        );

        //- Constructor for CAD files. Triangulation is done on the fly by OpenCASCADE
        triSurface
        (
            const fileName&,
            const scalar& linDeflection,
            const scalar& angDeflection,
            const bool forceTriangulation = false,
            const word nameRegionBy = "solid",
            const wordHashSet& ignoreRegions = wordHashSet()
        );

        //- Construct from Istream
        triSurface(Istream&);

        //- Construct from objectRegistry
        triSurface(const Time& d);

        //- Construct as copy
        triSurface(const triSurface&);


    //- Destructor
    virtual ~triSurface();

        void clearOut();

        void clearTopology();

        void clearPatchMeshAddr();


    // Member Functions

      // Access

            //- Sort faces according to region. Returns patch list
            //  and sets faceMap to index of labelledTri inside *this.
            surfacePatchList calcPatches(labelList& faceMap) const;

        const geometricSurfacePatchList& patches() const
        {
            return patches_;
        }

        geometricSurfacePatchList& patches()
        {
            return patches_;
        }

        //- Return const access to the faces
        inline const List<labelledTri>& surfFaces() const
        {
            return static_cast<const List<labelledTri>&>(*this);
        }

        //- Return edge-face addressing sorted (for edges with more than
        //  2 faces) according to the angle around the edge.
        //  Orientation is anticlockwise looking from
        //  edge.vec(localPoints())
        const labelListList& sortedEdgeFaces() const;

        //- If 2 face neighbours: label of face where ordering of edge
        //  is consistent with righthand walk.
        //  If 1 neighbour: label of only face.
        //  If >2 neighbours: undetermined.
        const labelList& edgeOwner() const;


        //- Face area vectors (normals)
        inline const vectorField& Sf() const
        {
            return ParentType::faceAreas();
        }

        //- access old points after stl has moved
        const pointField& points0() const;

        //- Face area magnitudes
        inline const scalarField& magSf() const
        {
            return ParentType::magFaceAreas();
        }

        //- Face centres
        inline const vectorField& Cf() const
        {
            return ParentType::faceCentres();
        }

        //- Calculate tri-surface curvature
        scalarField calcCurvature
        (
            const word tag,
            const label nSmoothIter
        );

      // Interoperability with other surface mesh classes

        //- Sort faces according to zoneIds
        //  Returns a surfZoneList and sets faceMap to index within faces()
        //  (i.e. map from original,unsorted to sorted)
        List<surfZone> sortedZones(labelList& faceMap) const;

        //- Create a list of faces from the triFaces
        void triFaceFaces(List<face>& plainFaceList) const;


      // Edit

        //- Move points
        virtual void movePoints(const pointField&);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar);

            //- Scale points. A non-positive factor is ignored
            virtual void scalePoints(const vector);

            //- Translate points. A non-positive factor is ignored
            virtual void translatePoints(const vector);

            //- rotate points.
            virtual void rotate(Pair<vector> n1n2);

            //- rotate points based on quaternion.
            virtual void rotate(quaternion R);

        //- Check/remove degenerate triangles
        void checkDegenerate(const bool verbose);

        //- Check/remove duplicate/degenerate triangles
        void checkTriangles(const bool verbose);

        //- Check triply (or more) connected edges.
        void checkEdges(const bool verbose);

        //- Remove non-valid triangles
        void cleanup(const bool verbose);

        //- Fill faceZone with currentZone for every face reachable
        //  from facei without crossing edge marked in borderEdge.
        //  Note: faceZone has to be sized nFaces before calling this fun.
        void markZone
        (
            const boolList& borderEdge,
            const label facei,
            const label currentZone,
            labelList& faceZone
        ) const;

        //- (size and) fills faceZone with zone of face. Zone is area
        //  reachable by edge crossing without crossing borderEdge
        //  (bool for every edge in surface). Returns number of zones.
        label markZones
        (
            const boolList& borderEdge,
            labelList& faceZone
        ) const;

        //- 'Create' sub mesh, including only faces for which
        //  boolList entry is true
        //  Sets: pointMap: from new to old localPoints
        //        faceMap: new to old faces
        void subsetMeshMap
        (
            const boolList& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface. Returns pointMap, faceMap from
        //  subsetMeshMap
        triSurface subsetMesh
        (
            const boolList& include,
            labelList& pointMap,
            labelList& faceMap,
            bool compactPatches = false
        ) const;


        //- Transfer stored faces to an Xfer container
        Xfer<List<labelledTri>> xferFaces();

        //- Transfer stored points to an Xfer container
        Xfer<List<point>> xferPoints();

        //- Alter contents by transferring (triangles, points) components.
        //  Patch information is small and therefore just copied.
        void transfer(triSurface& surf);

      // Write

        //- Write to Ostream in simple FOAM format
        void write(Ostream&) const;

        //- Generic write routine. Chooses writer based on extension.
        void write
        (
            const fileName&,
            const bool sortByRegion = false,
            const bool keepEmpty = false
        ) const;

        //- Write to database
        void write(const Time&) const;

        //- Write some statistics
        void writeStats(Ostream& os) const;


    // Member operators

        void operator=(const triSurface&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, triSurface&);
        friend Ostream& operator<<(Ostream&, const triSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
