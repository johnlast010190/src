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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::searchableSurfaces

Description
    Container for searchableSurfaces.

SourceFiles
    searchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaces_H
#define searchableSurfaces_H

#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "primitives/Pair/labelPair.H"
#include "sampledSetWriters/writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class searchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaces
:
    public PtrList<searchableSurface>
{
    // Private data

        //- Surface names
        wordList names_;

        //- Surface file names
        wordList fileNames_;

        //- Region names per surface
        List<wordList> regionNames_;

        //- Indices of all surfaces. Precalculated and stored.
        labelList allSurfaces_;

        //- Whether surface has been set as single region.
        boolList singleRegion_;

    // Private Member Functions

        //- Is edge on face
        static bool connected
        (
            const triSurface& s,
            const label edgeI,
            const pointIndexHit& hit
        );


        //- Disallow default bitwise copy construct
        searchableSurfaces(const searchableSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaces&);


public:

    ClassName("searchableSurfaces");


    // Constructors

        //- Construct with length specified. Fill later.
        explicit searchableSurfaces(const label);


        ////- Construct from list of dictionaries
        //searchableSurfaces(const IOobject&, const PtrList<dictionary>&);

        //- Construct from dictionary
        searchableSurfaces
        (
            const IOobject&,
            const dictionary&,
            const bool singleRegionName = false,
            const List<fileName>& alternativeSurfacePath = List<fileName>(0),
            const label maxNumberRegions = labelMax
        );


    // Member Functions

        const wordList& names() const
        {
            return names_;
        }

        wordList& names()
        {
            return names_;
        }

        const wordList& fileNames() const
        {
            return fileNames_;
        }
        wordList& fileNames()
        {
            return fileNames_;
        }

        const List<wordList>& regionNames() const
        {
            return regionNames_;
        }

        List<wordList>& regionNames()
        {
            return regionNames_;
        }

        const List<bool>& singleRegion() const
        {
            return singleRegion_;
        }
        List<bool>& singleRegion()
        {
            return singleRegion_;
        }


        ////- If adding surfaces 'by hand'
        //HashTable<labelPair>& regionNames()
        //{
        //    return regionNames_;
        //}
        ////- Get surface and region for a name
        //const labelPair& surfaceRegion(const word& globalRegion) const
        //{
        //    return regionNames_[globalRegion];
        //}

        //- Find index of surface. Return -1 if not found.
        label findSurfaceID(const word& name) const;

        //- Find index of surface from fileName. Return -1 if not found.
        label findFileSurfaceID(const word& name) const;

        label findSurfaceRegionID
        (
            const word& surfaceName,
            const word& regionName
        ) const;

        // Multiple point queries.

            //- Find any intersection. Return hit point information and
            //  surface number. If multiple surfaces hit the first surface
            //  is returned, not necessarily the nearest (to start).
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find all intersections in order from start to end. Returns for
            //  every hit the surface and the hit info.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                labelListList& surfaces,
                List<List<pointIndexHit>>&
            ) const;

            //Find intersections of edge nearest to both endpoints.
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& surface2,
                List<pointIndexHit>& hit2
            ) const;

            //- Find nearest. Return -1 (and a miss()) or surface and nearest
            //  point.
            void findNearest
            (
                const pointField&,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            void findNearest
            (
                const labelListList& regionIndices,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearestSurfaces,
                List<pointIndexHit>& nearestInfo
            ) const;

            //- Calculate bounding box
            boundBox bounds() const;

        // Checking

            //- Are all surfaces closed and manifold
            bool checkClosed(const bool report) const;

            //- Are all (triangulated) surfaces consistent normal orientation
            bool checkNormalOrientation(const bool report) const;

            //- Are all bounding boxes of similar size
            bool checkSizes(const scalar maxRatio, const bool report) const;

            //- Do surfaces self-intersect or intersect others
            bool checkIntersection
            (
                const scalar tol,
                const autoPtr<writer<scalar>>&,
                const bool report
            ) const;

            //- Check triangle quality
            bool checkQuality
            (
                const scalar minQuality,
                const bool report
            ) const;

            //- All topological checks. Return number of failed checks
            label checkTopology(const bool report) const;

            //- All geometric checks. Return number of failed checks
            label checkGeometry
            (
                const scalar maxRatio,
                const scalar tolerance,
                const autoPtr<writer<scalar>>& setWriter,
                const scalar minQuality,
                const bool report
            ) const;

            //- Write some stats
            void writeStats(const List<wordList>&, Ostream&) const;


    // Member Operators

        //- Return const and non-const reference to searchableSurface by index.
        using PtrList<searchableSurface>::operator[];

        //- Return const reference to searchableSurface by name.
        const searchableSurface& operator[](const word&) const;

        //- Return reference to searchableSurface by name.
        searchableSurface& operator[](const word&);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
