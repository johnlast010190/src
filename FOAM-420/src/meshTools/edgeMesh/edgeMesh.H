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
    (c) 2017 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::edgeMesh

Description
    Points connected by edges.

    Can be read from fileName based on extension. Uses ::New factory method
    to select the reader and transfer the result.

SourceFiles
    edgeMeshI.H
    edgeMesh.C
    edgeMeshIO.C
    edgeMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMesh_H
#define edgeMesh_H

#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/edge/edgeList.H"
#include "edgeMesh/edgeMeshFormats/edgeMeshFormatsCore.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/runTimeSelection/memberFunctions/memberFunctionSelectionTables.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "primitives/quaternion/quaternion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class edgeMesh;
Istream& operator>>(Istream& is, edgeMesh& em);
Ostream& operator<<(Ostream& os, const edgeMesh& em);


/*---------------------------------------------------------------------------*\
                           Class edgeMesh Declaration
\*---------------------------------------------------------------------------*/

class edgeMesh
:
    public fileFormats::edgeMeshFormatsCore
{
    // Private data

        //- Vertices of the edges
        pointField points_;

        //- The edges defining the boundary
        edgeList edges_;

        //- From point to edges
        mutable autoPtr<labelListList> pointEdgesPtr_;


    // Private Member Functions

        //- Calculate point-edge addressing (inverse of edges)
        void calcPointEdges() const;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        inline pointField& storedPoints();

        //- Non-const access to the edges
        inline edgeList& storedEdges();


public:

        //- Runtime type information
        TypeName("edgeMesh");


    // Static

        //- Can we read this file format?
        static bool canRead(const fileName& name, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format type?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        edgeMesh();

        //- Construct from components
        edgeMesh(const pointField& points, const edgeList& edges);

        //- Construct by transferring components (points, edges).
        edgeMesh
        (
            const Xfer<pointField>& pointLst,
            const Xfer<edgeList>& edgeLst
        );

        //- Construct as copy
        edgeMesh(const edgeMesh& em);

        //- Construct from file name (uses extension to determine type)
        edgeMesh(const fileName& name);

        //- Construct from file name (uses extension to determine type)
        edgeMesh(const fileName& name, const word& ext);

        //- Construct from Istream
        edgeMesh(Istream&);

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            edgeMesh,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<edgeMesh> New
        (
            const fileName& name,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<edgeMesh> New(const fileName& name);


    //- Destructor
    virtual ~edgeMesh();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            edgeMesh,
            write,
            fileExtension,
            (
                const fileName& name,
                const edgeMesh& mesh
            ),
            (name, mesh)
        );

        //- Write to file
        static void write(const fileName& name, const edgeMesh& mesh);


    // Member Functions

        //- Transfer the contents of the argument and annul the argument
        void transfer(edgeMesh&);

        //- Transfer contents to the Xfer container
        Xfer<edgeMesh> xfer();

    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& ext);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);

        //- Read in Foam format
        bool read(Istream&);

    // Access

        //- Return points
        inline const pointField& points() const;

        //- Return edges
        inline const edgeList& edges() const;

        //- Return edges
        inline const labelListList& pointEdges() const;

        //- Find connected regions. Set region number per edge.
        //  Returns number of regions.
        label regions(labelList& edgeRegion) const;


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Reset primitive data (points, edges)
        //  Note, optimized to avoid overwriting data (with Xfer::null)
        virtual void reset
        (
            const Xfer<pointField>& points,
            const Xfer<edgeList>& edges
        );

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar scaleFactor);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const vector scaleFactor);

        //- Translate points. A non-positive factor is ignored
        virtual void translatePoints(const vector translateVector);

        //- rotate points
        virtual void rotate(quaternion R);

        //- rotate points based on quaternion
        virtual void rotate(Pair<vector> n1n2);

        //- Perform transforms based on dictionary input
        void doTransforms(const dictionary& dict);

        //- Geometric merge points (points within mergeDist) prior to
        //  automatically calling mergeEdges().
        virtual void mergePoints(const scalar mergeDist);

        //- Merge duplicate edges and eliminate unused points.
        virtual void mergeEdges();

        void resetEdgeList(edgeList &newEdges) { edges_.transfer(newEdges); }

        //- Find points used in subset
        static void subsetEdgeMeshMap
        (
            const edgeMesh& eMesh,
            const boolList& include,
            const label nIncluded,
            labelList& newToOldPoints,
            labelList& oldToNewPoints,
            labelList& newToOldEdges
        );

        //- Construct subsetted edgeMesh
        static edgeMesh subsetEdgeMesh
        (
            const edgeMesh& eMesh,
            const labelList& newToOldPoints,
            const labelList& oldToNewPoints,
            const labelList& newToOldEdges
        );

        //- Subset given marked edges
        static edgeMesh subsetEdgeMesh
        (
            const edgeMesh& eMesh,
            const boolList& include,
            labelList& newToOldPoints,
            labelList& newToOldEdges
        );

        //- Subset given marked edges
        static edgeMesh subsetEdgeMesh
        (
            const edgeMesh& eMesh,
            const labelList& newToOldEdges,
            labelList& newToOldPoints
        );

        //- Find edge otherE in allEdges.
        static label findEdge
        (
            const List<edge>& allEdges,
            const labelListList& allPointEdges,
            const edge& otherE
        );

        //- Merge edgeMesh (subEdges, subPoints) into allEdges, allPoints.
        static void merge
        (
            const scalar mergeDist,
            const List<edge>& subEdges,
            const pointField& subPoints,

            List<edge>& allEdges,
            pointField& allPoints,

            labelList& edgeConstructMap,
            labelList& pointConstructMap
        );

        //- Subset the part of edgeMesh that is overlapping bounds.
        static edgeMesh overlappingEdgeMesh
        (
            const edgeMesh&,
            const List<treeBoundBox>&,
            labelList& subPointMap,
            labelList& subEdgeMap
        );

        //- Does any part of edge overlap bb.
        static bool overlaps
        (
            const List<treeBoundBox>& bb,
            const point& p0,
            const point& p1
        );

        //- reconstruct all edges meshes
        void reconstruct(const scalar& mergeDist);

        //- Trim edges not contained in bounding box
        void trim
        (
            const List<treeBoundBox>& bbs
        );

        void distribute
        (
            const List<treeBoundBox>& bbs,
            const scalar& mergeDist
        );

    // Write

        virtual void writeStats(Ostream&) const;

        //- Generic write routine. Chooses writer based on extension.
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }

    // Member Operators

        inline void operator=(const edgeMesh& rhs);

        // Ostream Operator

            friend Ostream& operator<<(Ostream& os, const edgeMesh& em);
            friend Istream& operator>>(Istream& is, edgeMesh& em);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeMesh/edgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
