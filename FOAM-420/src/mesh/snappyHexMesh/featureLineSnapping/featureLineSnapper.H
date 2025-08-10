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
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::featureLineSnapper

Description
    Helper class which ensures specified geometry feature lines appear in the
    snapped surface mesh.

SourceFiles
    featureLineSnapper.C

\*---------------------------------------------------------------------------*/

#ifndef featureLineSnapper_H
#define featureLineSnapper_H

#include <utility>

#include "primitives/Vector/vector/vector.H"
#include "containers/Lists/List/ListDecl.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "memory/autoPtr/autoPtr.H"
#include "indexedOctree/treeDataEdge.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "featureLineSnapping/featureLinePrep.H"
#include "primitives/Pair/labelPair.H"
#include "patchDist/wallPoint/wallPoint.H"

namespace Foam
{
    class fvMesh;
    class mapPolyMesh;
    class polyTopoChange;
    class meshRefinement;
    class dictionary;
    class featureEdgeMesh;

    typedef featureLinePrep::featureLine featureLine; //just for convenience


    //- Find a unique edge direction.
    struct uniqueEdgeDir
    {
        public:

        void operator()(vector& x, const vector& y) const
        {
            if (x != vector(GREAT, GREAT, GREAT))
            {
                return;
            }
            else if (y != vector(GREAT, GREAT, GREAT))
            {
                x = y;
            }
        }
    };

    // To combineReduce a labelList. Filters out duplicates.
    struct uniqueEqOp
    {
       void operator()(labelList& x, const labelList& y) const
       {
          if (x.empty())
          {
             if (y.size())
             {
                x = y;
             }
          }
          else
          {
             forAll(y, yi)
             {
                if (findIndex(x, y[yi]) == -1)
                {
                   label sz = x.size();
                   x.setSize(sz+1);
                   x[sz] = y[yi];
                }
             }
          }
       }
    };

    //- Combine operator class for equalizing displacements.
    struct minMagEqOp
    {
        void operator()(vector& x, const vector& y) const
        {
            if (magSqr(y) < magSqr(x)) {x = y;}
        }
    };

    struct maxMagEqOp
    {
        void operator()(vector& x, const vector& y) const
        {
            if (magSqr(y) > magSqr(x)) {x = y;}
        }
    };

    ///Creates a displacement field for snapping points to feature lines
    /**
     * The supplied featureLinePrep object contains feature lines extracted
     * from the geometry.  This class implements an algorithm to build their
     * "shadows" in the corresponding mesh boundary patches, thereby
     * constructing a patch displacement field which can be applied to the
     * patch as part of the larger feature snapping algorithm in autoSnapDriver.
     */
    class featureLineSnapper
    {
        //nested classes

        /// private class for implementing the node-manifold-to-feature
        //line projection
        struct nodeMfdProjResult
        {
            nodeMfdProjResult() :
            index(-1),
            distSqr(GREAT),
            hit(pointIndexHit())
            {}
            nodeMfdProjResult
            (
                label index,
                scalar distSqr,
                const pointIndexHit &hit
            ) :
            index(index),
            distSqr(distSqr),
            hit(hit)
            {}

            label index; //index in node manifold
            scalar distSqr; //squared distance from feature line
            pointIndexHit hit; //hit object returned by tree
        };

        ///Sister class to nodeMfdProjResult, defining a partial order based
        ///on snapping distance
        struct nMPRCmp
        {
            bool operator()
            (
                const nodeMfdProjResult &a,
                const nodeMfdProjResult &b
            )
            {
                return a.distSqr < b.distSqr;
            }
        };

        //- Reference data from outside
        const featureLinePrep &flPrep_;

        const dictionary &motionDict_;

        const labelList &adaptPatchIDs_;

        indirectPrimitivePatch &pp_;

        meshRefinement &meshRefiner_;

        List<labelPair>& baffles_;

        const bool dualMesh_;

        const bool allowQuadSplits_;

        //- Internal data

        //- Displacements of patch points (mostly zero, of course)
        pointField pointDisp_;

        //- The chain data themselves - lists of (local) node numbers
        List<labelList> shadowChains_;

        //- The raw chain displacement data
        List<List<vector>> shadowChainDisp_;

        //- Whether a surface boundary feature
        List<bool> shadowChainSurfBound_;

        //- Array to flag faces that have already been split by some top mod
        //- op and mustn't therefore be split again
        PackedList<1> faceIsSplit;

        scalar maxShrinkFactor_;  //amount by which fl snapping is allowed to
                  //shrink edges and min face angles before it
                  //is cancelled

        scalar minFaceAngle_;  //check min face angle for splitting

        scalar minDihedralForSplits_; //do we allow the chain to follow face
            //diagonals?

        //- Debug configuration
        bool verbose_; //trigger lots of messages from the chaining algorithm
        bool dump_; //dump vtk files of the features lines and shadow chains

        mutable label nNoValidEndPts_; //count of how many feature lines were
                                   //too far from the patch to have a shadow
                           //chain constructed safely
        mutable label nFailedConstr_; //count of how many shadow chain
                          //constructions were abandoned because a
                          //new valid node couldn't be found

    public:

        featureLineSnapper
        (
            const featureLinePrep &flPrep,
            const dictionary &snapDict,
            const dictionary &motionDict,
            const labelList &adaptPatchIDs,
            indirectPrimitivePatch &pp,
            meshRefinement &meshRefiner,
            List<labelPair>& baffles,
            bool dualMesh,
            bool allowQuadSplits
        );


        //- Attempts to build patch-point chains that shadow feature lines
        void constructShadowChains
        (
            const PackedList<1>& isZonedFace
        );

        //Updates after topology changes
        void updateMesh(const mapPolyMesh& map);

        //- Renumbers chains and displacement field if patch numbering is
        //changed for any reason
        void renumber(const labelList &localPointMap);

        //- performs snapping of the shadow chain points to the corresponding
        //- feature line points
        void snapShadowChainsToLines
        (
            pointField &newPoints,
            boolList &stationaryPoints
        ) const;

        //- simple function building up a flag list indicating which nodes
        //are in chains
        autoPtr<PackedList<1>> flagChainNodes() const;

        List<labelList> getMeshNumberedShadowChains() const;

        const List<labelList>& getShadowChains() const
        {
            return shadowChains_;
        }

        pointField& pointDisp(){ return pointDisp_; }

    private:


        //- Return list of neighbouring points
        labelList neigbouringPoints
        (
            const label curPt,
            const bool allowSplits,
            const bool includeCurrentPt
        ) const;

        //- Find patch points corresponding to start and end points
        //- of feature line
        void findEndPts
        (
            const treeBoundBox& bb,
            const featureLine& featurel,
            const List<scalar>& flParams,
            const indexedOctree<treeDataEdge>& tree,
            const List<scalar>& lengthAve,
            label& minStartPt,
            label& minEndPt
        );

        //- loop through the feature lines constructing their "shadows" on the
        //- patch surface
        void buildShadowChains();

        //- simplify chains where they partially coincide
        label filterChains();

        void chainDispToPatchDisp();

        //-remove two-point singularities (nodes connected to only two edges)
        //-from the patch in order to improve the quality of chain-snapped mesh
        void removeTopoSingularities(fvMesh &mesh,polyTopoChange &meshMod);

        //- loop through getting rid of unwanted chain topologies
        void cleanUpChains();

        //-zero displacements not corresponding to a shadow chain node
        void cleanDisplacementField();

        //- Pull back any displacements that are greater than
        //the local patch edge length
        void filterDisplacementField();

        //- Apply a little bit of smoothing to non-chain nodes
        void smoothDisplacementField();

        //- Shouldn't really be in this class - just does a simple calculation
        //- based on the points defining a polygon
        static scalar getMinFaceAngle(const List<vector> &points);

        //- modify the mesh by splitting faces where the chain goes
        //across a diagonal
        void addChainDiagonalEdges
        (
            const PackedList<1>& isZonedFace
        );

        //- find the next point in the shadow chain for the current feature line
        //- return 0 if a valid next point other than the end point was found
        //- 1 if the end point was found
        //- 2 if the search failed
        label findNext
        (
            label &curPt,
            scalar &curFtrDist,
            vector &curTargetPt,
            const indexedOctree<treeDataEdge> &ftrLineTree,
            const featureLine &ftrLine,
            const List<scalar> &flParams,
            const bool procPt,
            const label endPt,
            const List<scalar> &lengthAve,
            vectorField& allDisp,
            bool reverse
        ) const;

        bool chainIsValid
        (
            DynamicList<label> &chain,
            DynamicList<vector> &chainDisp,
            const List<scalar> &lengthAve
         ) const;

        //- check whether the two points define an edge in the surface patch
        //- (assumes local point numbering)
        bool isEdge(label pt1, label pt2) const;

        //- set up the face split defined by the two points (local numbering)
        //- which are assumed to lie diagonally opposite on a face.  Returns
        //- false if in fact they don't.
        bool getFaceSplit
        (
            label pt1,
            label pt2,
            label &origFace,
            face &origNew,
            face &newNew
        ) const;

    };

}
#endif
