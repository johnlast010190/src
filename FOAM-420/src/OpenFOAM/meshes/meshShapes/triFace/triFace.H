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

Class
    Foam::triFace

Description
    A triangular face using a FixedList of labels corresponding to mesh
    vertices.

See also
    Foam::face, Foam::triangle

SourceFiles
    triFaceI.H
    triFaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef triFace_H
#define triFace_H

#include "containers/Lists/FixedList/FixedList.H"
#include "meshes/meshShapes/edge/edgeList.H"
#include "meshes/primitiveShapes/objectHit/pointHit.H"
#include "meshes/primitiveShapes/triangle/intersection.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/primitiveShapes/triangle/triPointRef.H"
#include "containers/Lists/ListListOps/ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class face;
class triFace;

inline bool operator==(const triFace& a, const triFace& b);
inline bool operator!=(const triFace& a, const triFace& b);


/*---------------------------------------------------------------------------*\
                          Class triFace Declaration
\*---------------------------------------------------------------------------*/

class triFace
:
    public FixedList<label, 3>
{

public:

    // Constructors

        //- Construct null with invalid point labels (-1)
        inline triFace();

        //- Construct from three point labels
        inline triFace
        (
            const label a,
            const label b,
            const label c
        );

        //- Construct from a list of 3 labels.
        explicit inline triFace(const labelUList& lst);

        //- Construct from an initializer list of 3 labels
        explicit inline triFace(std::initializer_list<label> lst);

        //- Construct from Istream
        inline triFace(Istream& is);


    // Member Functions

        //- 'Collapse' face by marking duplicate point labels.
        //  Duplicates point labels are marked with '-1'
        //  (the lower vertex is retained).
        //  Return the collapsed size.
        inline label collapse();

        //- Flip the face in-place.
        //  The starting points of the original and flipped face are identical.
        inline void flip();

        //- Return the points corresponding to this face
        inline pointField points(const UList<point>& points) const;

        //- Return triangle as a face
        inline face triFaceFace() const;

        //- Return the triangle
        inline triPointRef tri(const UList<point>& points) const;

        //- Return centre (centroid)
        inline point centre(const UList<point>& points) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average(const UList<point>& unused, const Field<Type>& fld) const;

        //- Magnitude of face area
        inline scalar mag(const UList<point>& points) const;

        //- The area normal - with magnitude equal to area of face
        inline vector areaNormal(const UList<point>& points) const;

        //- The unit normal
        inline vector unitNormal(const UList<point>& points) const;

        //- Legacy name for areaNormal()
        //  \deprecated(v3.2) Deprecated for new use
        FOAM_DEPRECATED_FOR(v3.2, "areaNormal() or unitNormal()")
        vector normal(const UList<point>& points) const
        {
            return areaNormal(points); // Legacy definition
        }

        //- Number of triangles after splitting
        inline label nTriangles() const;

        //- Return face with reverse direction
        //  The starting points of the original and reverse face are identical.
        inline triFace reverseFace() const;

        //- Return true if the global point label is found in face.
        bool found(const label globalIndex) const;

        //- Which local index (0,1,2) on face given a global index.
        //  returns -1 if not found
        label which(const label globalIndex) const;

        //- Return swept-volume from old-points to new-points
        inline scalar sweptVol
        (
            const UList<point>& opts,
            const UList<point>& npts
        ) const;

        //- Return the inertia tensor, with optional reference
        //  point and density specification
        inline tensor inertia
        (
            const UList<point>& points,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return point intersection with a ray starting at p, in direction q.
        inline pointHit ray
        (
            const point& p,
            const vector& q,
            const UList<point>& points,
            const intersection::algorithm = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        inline pointHit intersection
        (
            const point& p,
            const vector& q,
            const UList<point>& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        inline pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const UList<point>& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        inline pointHit nearestPoint
        (
            const point& p,
            const UList<point>& points
        ) const;


        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        inline pointHit nearestPointClassify
        (
            const point& p,
            const UList<point>& points,
            label& nearType,
            label& nearLabel
        ) const;

        //- based on the point coordinated of the triangle,
        // and a point inside, find the weights which the point can
        // be reconstructed based on the points
        // p ~ pRec = w0*p0 + w1*p1 + w2*p2;
        inline vector findTriangleWeights
        (
            const point& p,
            const pointField& points
        ) const;

        //- Return number of edges
        inline label nEdges() const;

        //- Return edges in face point ordering,
        //  i.e. edges()[0] is edge between [0] and [1]
        inline edgeList edges() const;

        //- Walk the edges in point order. Use with a range-for loop.
        face_edge_walker<triFace> walkEdges() const {
            return face_edge_walker<triFace>{*this};
        }

        //- Return n-th face edge
        inline edge faceEdge(const label n) const;

        //- Return the edge direction on the face
        //  Returns:
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        //  -  0: edge not found on the face
        inline int edgeDirection(const edge& e) const;

        //- Compare triFaces
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same face, but different orientation
        static inline int compare(const triFace& a, const triFace& b);

};


template<>
constexpr bool contiguous<triFace>()  {return true;}

inline bool operator==(const triFace& a, const triFace& b);
inline bool operator!=(const triFace& a, const triFace& b);


//- Hash specialization for hashing triFace - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<triFace>::operator()(const triFace& t, unsigned seed) const
{
    // Fortunately we don't need this very often
    const uLabel t0(t[0]);
    const uLabel t1(t[1]);
    const uLabel t2(t[2]);

    const uLabel val = (t0*t1*t2 + t0+t1+t2);

    return Hash<uLabel>()(val, seed);
}


//- Hash specialization for hashing triFace - a commutative hash value.
template<>
inline unsigned Hash<triFace>::operator()(const triFace& t) const
{
    return Hash<triFace>::operator()(t, 0);
}


//- Hash specialization to offset faces in ListListOps::combineOffset
template<>
class offsetOp<triFace>
{

public:

    inline triFace operator()
    (
        const triFace& x,
        const label offset
    ) const
    {
        triFace result(x);

        forAll(x, xI)
        {
            result[xI] = x[xI] + offset;
        }
        return result;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshes/meshShapes/triFace/triFaceI.H"

#ifdef NoRepository
    #include "meshes/meshShapes/triFace/triFaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
