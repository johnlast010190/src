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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::triSurfaceMesh

Description
    IOoject and searching on triSurface

    Note: when constructing from dictionary has optional parameters:
        - scale     : scaling factor.
        - tolerance : relative tolerance for doing intersections
                      (see triangle::intersection)
        - minQuality: discard triangles with low quality when getting normal

SourceFiles
    triSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMesh_H
#define triSurfaceMesh_H

#include "triSurface/triSurfaceSearch/triSurfaceRegionSearch.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "db/objectRegistry/objectRegistry.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataTriSurface.H"
#include "indexedOctree/treeDataPrimitivePatch.H"
#include "indexedOctree/treeDataEdge.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class triSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMesh
:
    public searchableSurface,
    public objectRegistry,      // so we can store fields
    public triSurface,
    public triSurfaceRegionSearch
{
protected:

    // Private member data

        //- Supplied fileName override
        fileName fName_;

        //- Optional min triangle quality. Triangles below this get
        //  ignored for normal calculation
        scalar minQuality_;

        //- Search tree for boundary edges.
        mutable autoPtr<indexedOctree<treeDataEdge>> edgeTree_;

        //- Names of regions
        mutable wordList regions_;

        //- Is surface closed
        mutable label surfaceClosed_;

        //- If surface is closed, what is type of outside points
        mutable volumeType outsideVolType_;


    // Private Member Functions

        //- Return fileName to load IOobject from
        static fileName checkFile(const regIOobject& io, const bool isGlobal);

        //- Return fileName to load IOobject from. Optional override of fileName
        static fileName checkFile
        (
            const regIOobject&,
            const dictionary&,
            const bool isGlobal
        );

        //- Helper function for isSurfaceClosed
        static bool addFaceToEdge
        (
            const edge&,
            EdgeMap<label>&
        );

        //- Check whether surface is closed without calculating any permanent
        //  addressing.
        bool isSurfaceClosed() const;

        //- Steps to next intersection. Adds smallVec and starts tracking
        //  from there.
        static void getNextIntersections
        (
            const indexedOctree<treeDataTriSurface>& octree,
            const point& start,
            const point& end,
            const vector& smallVec,
            DynamicList<pointIndexHit, 1, 1>& hits
        );

        //- Disallow default bitwise copy construct
        triSurfaceMesh(const triSurfaceMesh&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMesh&);


public:

    //- Runtime type information
    TypeName("triSurfaceMesh");

    //- Return the mesh sub-directory name (usually "triSurface")
    static word meshSubDir;


    // Constructors

        //- Construct from triSurface
        triSurfaceMesh(const IOobject&, const triSurface&);

        //- Construct read
        triSurfaceMesh(const IOobject& io);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        triSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


        // Special constructors for use by distributedTriSurface. File search
        // method supplied:

            enum readAction
            {
                localOnly,      // load from (processor-)local file only
                localOrGlobal,  // load from (processor-)local or global file
                masterOnly      // as localOrGlobal but only load on master
            };

            triSurfaceMesh(const IOobject& io, const readAction r);

            triSurfaceMesh
            (
                const IOobject& io,
                const dictionary& dict,
                const readAction r
            );


    //- Destructor
    virtual ~triSurfaceMesh();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Move points
        virtual void movePoints(const pointField&);

        //- Demand driven construction of octree for boundary edges
        const indexedOctree<treeDataEdge>& edgeTree() const;


        // searchableSurface implementation

            virtual const wordList& regions() const;

            //- Whether supports volume type below. I.e. whether is closed.
            virtual bool hasVolumeType() const;

            //- Range of local indices that can be returned.
            virtual label size() const
            {
                return triSurface::size();
            }

            //- Get representative set of element coordinates
            //  Usually the element centres (should be of length size()).
            virtual tmp<pointField> coordinates() const;

            //- Get bounding spheres (centre and radius squared). Any point
            //  on surface is guaranteed to be inside.
            virtual void boundingSpheres
            (
                pointField& centres,
                scalarField& radiusSqr
            ) const;

            //- Get the points that define the surface.
            virtual tmp<pointField> points() const;

            // Does any part of the surface overlap the supplied bound box?
            virtual bool overlaps(const boundBox& bb) const;

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- Find all elements intersecting sphere
            virtual void findSphere
            (
                const pointField& centres,
                const scalarField& radiusSqr,
                List<labelList>& hitIndexes,
                labelList& nOverlapChecks
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region,
                const bool threaded = false
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&,
                const bool threaded = false
            ) const;


        // Other

            //- WIP. Store element-wise field.
            virtual void setField(const labelList& values);

            //- WIP. From a set of hits (points and
            //  indices) get the specified field. Misses do not get set.
            virtual void getField
            (
                const List<pointIndexHit>&,
                labelList&,
                const bool threaded = false
            ) const;

            //- WIP. Clear out curvature field if registered.
            virtual void clearCurvaturePointField();

            //- WIP. Store point-wise curvature field.
            virtual void setCurvaturePointField(const scalarField& curvature);

            //- WIP. From a set of hits (points and
            //  indices) get the specified curvature field. Misses do not get set.
            virtual void getCurvatureField(const List<pointIndexHit>&, scalarField&) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
            }

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool valid
            ) const;

            //- Is object global
            virtual bool global() const
            {
                return true;
            }

            //- Return complete path + object name if the file exists
            //  either in the case/processor or case otherwise null
            virtual fileName filePath() const
            {
                return searchableSurface::globalFilePath(type());
            }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<triSurfaceMesh>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
