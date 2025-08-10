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
    (c) 2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::face

Description
    A face is a list of labels corresponding to mesh vertices.

See also
    Foam::triFace

SourceFiles
    faceI.H
    face.C
    faceIntersection.C
    faceContactSphere.C
    faceAreaInContact.C
    faceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef face_H
#define face_H

#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/meshShapes/edge/edgeList.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "meshes/meshShapes/face/faceListFwd.H"
#include "meshes/primitiveShapes/triangle/intersection.H"
#include "meshes/primitiveShapes/objectHit/pointHit.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "faceEdgeIterator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class face;
class triFace;

inline Istream& operator>>(Istream& is, face& f);

/*---------------------------------------------------------------------------*\
                           Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public labelList
{
public:

    //- Return types for classify
    enum proxType
    {
        NONE,
        POINT,  // Close to point
        EDGE    // Close to edge
    };

    // Static data members

        static const char* const typeName;


    // Static Functions

        //- Return centre point given face points
        template<class PointField>
        static vector faceCentre(const PointField& ps);

        //- Return vector area given face points
        template<class PointField>
        static vector faceArea(const PointField& ps);


    // Constructors

        //- Construct null
        inline face();

        //- Construct given size, with invalid point labels (-1)
        explicit inline face(const label sz);

        //- Construct from list of labels
        explicit inline face(const labelUList& lst);

        //- Construct from an initializer list of labels
        explicit inline face(std::initializer_list<label> lst);

        //- Construct by transferring the parameter contents
        explicit inline face(const Xfer<labelList>& lst);

        //- Move constructor.
        face(face&& f) = default;
        face(const face& f) = default;
        face& operator=(const face& a) = default;
        face& operator=(face&& a) = default;

        //- Copy construct from triFace
        face(const triFace& f);

        //- Construct from Istream
        inline face(Istream& is);


    // Member Functions

        //- Collapse face by removing duplicate point labels
        //  return the collapsed size
        label collapse();

        //- Flip the face in-place.
        //  The starting points of the original and flipped face are identical.
        void flip();

        //- Return the points corresponding to this face
        inline pointField points(const UList<point>& points) const;

        //- Centre point of face
        point centre(const UList<point>& points) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average
        (
            const UList<point>& meshPoints,
            const Field<Type>& fld
        ) const;

        //- The area normal - with magnitude equal to area of face
        vector areaNormal(const UList<point>& p) const;

        //- The unit normal
        inline vector unitNormal(const UList<point>& p) const;

        //- Legacy name for areaNormal()
        //  \deprecated(v3.2) Deprecated for new use
        FOAM_DEPRECATED_FOR(v3.2, "areaNormal() or unitNormal()")
        vector normal(const UList<point>& p) const
        {
            return areaNormal(p); // Legacy definition
        }

        //- Magnitude of face area
        inline scalar mag(const UList<point>& p) const;

        //- Return face with reverse direction
        //  The starting points of the original and reverse face are identical.
        face reverseFace() const;

        //- Navigation through face vertices

            //- Return true if the global point label is found in face.
            inline bool found(const label globalIndex) const;

            //- Which local vertex on face given a global index.
            //  returns -1 if not found
            label which(const label globalIndex) const;

            //- Next vertex on face
            inline label nextLabel(const label i) const;

            //- Previous vertex on face
            inline label prevLabel(const label i) const;


        //- Return the volume swept out by the face when its points move
        scalar sweptVol
        (
            const UList<point>& oldPoints,
            const UList<point>& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //  point and density specification
        tensor inertia
        (
            const UList<point>& p,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return potential intersection with face with a ray starting
        //  at p, direction n (does not need to be normalized)
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. Face-centre is calculated from point average.
        //  For a hit, the distance is signed.  Positive number
        //  represents the point in front of triangle
        //  In case of miss the point is the nearest point on the face
        //  and the distance is the distance between the intersection point
        //  and the original point.
        //  The half-ray or full-ray intersection and the contact
        //  sphere adjustment of the projection vector is set by the
        //  intersection parameters
        pointHit ray
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. See triangle::intersection for details.
        pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const UList<point>& meshPoints,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        pointHit nearestPoint
        (
            const point& p,
            const UList<point>& meshPoints
        ) const;

        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        pointHit nearestPointClassify
        (
            const point& p,
            const UList<point>& meshPoints,
            label& nearType,
            label& nearLabel
        ) const;

        //- Return contact sphere diameter
        scalar contactSphereDiameter
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints
        ) const;

        //- Return area in contact, given the displacement in vertices
        scalar areaInContact
        (
            const UList<point>& meshPoints,
            const scalarField& v
        ) const;

        //- Return number of edges
        inline label nEdges() const;

        //- Walk the edges in point order. Use with a range-for loop.
        face_edge_walker<face> walkEdges() const {
            return face_edge_walker<face>{*this};
        }

        //- Return edges in face point ordering,
        //  i.e. edges()[0] is edge between [0] and [1]
        // Returns the result in a newly-allocated list, which is useful if you need to create
        // and store a copy of this list somewhere.
        //
        // If you just need to iterate the edges, you can save time by doing:
        //
        // for (const auto& edge : thisObject.walkEdges()) {
        //      // do something with `edge`
        // }
        //
        //
        // If you need to perform random access to edges, calling `faceEdge(n)` to access the one edge you want each
        // time will be significantly faster than `edges()[n]` (which allocates and then discards the entire list of
        // edges!).
        edgeList edges() const;

        //- Return n-th face edge
        __attribute__((always_inline))
        edge faceEdge(const label n) const;

        //- Return the edge direction on the face
        //  Returns:
        //  -  0: edge not found on the face
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        int edgeDirection(const edge& e) const;

        //- Find the longest edge on a face.
        label longestEdge(const UList<point>& pts) const;

        //- Size of the face's triangulation
        inline label nTriangles() const;

        //- Compare faces
        //   0: different
        //  +1: identical
        //  -1: same face, but different orientation
        static int compare(const face& a, const face& b);

        //- Return true if the faces have the same vertices
        static bool sameVertices(const face& a, const face& b);


    // Istream Operator

        friend Istream& operator>>(Istream& is, face& f);
};


//- Hash specialization to offset faces in ListListOps::combineOffset
template<>
class offsetOp<face>
{

public:

    inline face operator()
    (
        const face& x,
        const label offset
    ) const
    {
        face result(x.size());

        forAll(x, xI)
        {
            result[xI] = x[xI] + offset;
        }
        return result;
    }
};


// Global operators
inline bool operator==(const face& a, const face& b);
inline bool operator!=(const face& a, const face& b);


// Global functions

//- Find the longest edge on a face. Face point labels index into pts.
//  \deprecated use class method instead (APR-2017)
label longestEdge(const face& f, const UList<point>& pts);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshes/meshShapes/face/faceI.H"

#ifdef NoRepository
    #include "meshes/meshShapes/face/faceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
