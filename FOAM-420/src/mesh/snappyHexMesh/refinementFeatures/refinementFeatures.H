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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::refinementFeatures

Description
    Encapsulates queries for features.

SourceFiles
    refinementFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef refinementFeatures_H
#define refinementFeatures_H

#include "edgeMesh/extendedEdgeMesh/extendedFeatureEdgeMesh/extendedFeatureEdgeMesh.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataEdge.H"
#include "indexedOctree/treeDataPoint.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "fvMesh/fvMesh.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "refinementFeatures/manifoldFeatures.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class refinementFeatures Declaration
\*---------------------------------------------------------------------------*/

class refinementFeatures
:
    public PtrList<extendedFeatureEdgeMesh>
{
    // Private data

        // whether to check for refinement only features
        bool checkRefinementOnly_;

        //- Per shell if refinement feature only
        List<bool> refinementOnly_;

        //- Edge meshes subsetted into manifold feature edges
        manifoldFeatures mFeatures_;

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;

        //- Edge
        mutable PtrList<indexedOctree<treeDataEdge>> edgeTrees_;

        //- Features points
        mutable PtrList<indexedOctree<treeDataPoint>> pointTrees_;

        //- Region edge trees (demand driven)
        mutable autoPtr<PtrList<indexedOctree<treeDataEdge>>>
            regionEdgeTrees_;

        //- geometric angle along feature edge
        const scalar minCos_;

    // Private Member Functions

        //- Read set of feature edge meshes
        void read(const objectRegistry&, const PtrList<dictionary>&);

        //- Find shell level higher than ptLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label featI,
            labelList& maxLevel
        );


protected:

        const PtrList<indexedOctree<treeDataEdge>>& edgeTrees() const;

        const PtrList<indexedOctree<treeDataPoint>>& pointTrees() const;

        const PtrList<indexedOctree<treeDataEdge>>& regionEdgeTrees() const;

public:

    // Constructors

        //- Construct from description and refinemenet surfaces
        refinementFeatures
        (
            const fvMesh& mesh,
            const PtrList<dictionary>& featDicts,
            const refinementSurfaces& refineSurfaces,
            const refinementParameters& refineParams,
            const scalar& mergeDist
        );

        //- Construct from description
        refinementFeatures
        (
            const objectRegistry& io,
            const PtrList<dictionary>& featDicts
        );

        //- Construct from description and do geometric analysis to determine
        //  feature points
        refinementFeatures
        (
            const objectRegistry& io,
            const PtrList<dictionary>& featDicts,
            const scalar minCos
        );

    // Member Functions

        // Access

            const manifoldFeatures& manFeatures() const
            {
                return mFeatures_;
            }

            manifoldFeatures& manFeatures()
            {
                return mFeatures_;
            }

            //- Per edgeMesh whether only a refinement feature
            const List<bool>& refinementOnly() const
            {
                return refinementOnly_;
            }

            //- Per featureEdgeMesh the list of level
            const labelListList& levels() const
            {
                return levels_;
            }

            //- Per featureEdgeMesh the list of ranges
            const List<scalarField>& distances() const
            {
                return distances_;
            }

            void clearTrees();

            //- Set whether to check for refinement only surfaces
            void setCheckRefinementOnly(bool check)
            {
                checkRefinementOnly_ = check;
            }

            const bool& checkRefinementOnly() const
            {
                return checkRefinementOnly_;
            }

        // Helper

            //trim feature lines to mesh bounding box
            void trim(const polyMesh& mesh);


         // Query

            //- Highest distance of all features
            scalar maxDistance() const;

            //- Find nearest point on nearest feature edge. Sets
            //  - nearFeature: index of feature mesh
            //  - nearInfo   : location on feature edge and edge index
            //                 (note: not feature edge index but index into
            //                  edges() directly)
            //  - nearNormal : local feature edge normal
            void findNearestEdge
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo,
                vectorField& nearNormal
            ) const;

            //- Find nearest point on nearest region edge. Sets
            //  - nearFeature: index of feature mesh
            //  - nearInfo   : location on feature edge and edge index
            //                 (note: not feature edge index but index into
            //                  edges() directly)
            //  - nearNormal : local feature edge normal
            void findNearestRegionEdge
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo,
                vectorField& nearNormal
            ) const;

            //- Find nearest feature point. Sets
            //  - nearFeature: index of feature mesh
            //  - nearInfo   : location on feature point and point index.
            //                 (note: not index into shapes().pointLabels() but
            //                  index into points() directly)
            void findNearestPoint
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo
            ) const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            );

        // Other

            //- Check sizes - return true if error and stream to os
            bool checkSizes
            (
                const scalar maxRatio,
                const boundBox& meshBb,
                const bool report,
                Ostream& os
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
