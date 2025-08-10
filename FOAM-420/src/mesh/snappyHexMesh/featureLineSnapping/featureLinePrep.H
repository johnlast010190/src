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
    Foam::featureLinePrep

Description
    Helper class which ensures specified geometry feature lines appear in the
    snapped surface mesh.

SourceFiles
    featureLinePrep.C

\*---------------------------------------------------------------------------*/

#ifndef featureLinePrep_H
#define featureLinePrep_H

#include <utility>

#include "primitives/Vector/vector/vector.H"
#include "containers/Lists/List/ListDecl.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "memory/autoPtr/autoPtr.H"
#include "indexedOctree/treeDataEdge.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/Tuple2/Tuple2.H"
#include "edgeMesh/edgeMesh.H"
#include "refinementFeatures/refinementFeatures.H"

namespace Foam
{
    class fvMesh;
    class mapPolyMesh;
    class polyTopoChange;
    class meshRefinement;
    class dictionary;
    class featureEdgeMesh;
    class edgeMesh;
    class treeDataFace;
    class triSurfaceMesh;


    ///Constructs feature line objects from various sources.  This class is
    ///used as input to featureLineSnapper.
    class featureLinePrep
    {
        //nested classes

        public:
        /// List wrapper that gives easy access to coordinate information
        //  from various potential sources.
        class featureLine : public labelList
        {
            public:

            featureLine() : points_(0), minDihedral_(true) {}
            explicit featureLine
            (
                const pointField &points,
                scalar minDihedral
             )
             : points_(&points), minDihedral_(minDihedral)
             {}
            const vector & point(label j) const
            {
                return (*points_)[(*this)[j]];
            }
            const pointField &points() const
            {
                return *points_;
            }

            const edgeList &edges() const;

            scalar minDihedral() const
            {
                return minDihedral_;
            }

            private:

            //this is a reference, not ownership (just in case that's
            //not obvious)
            const pointField *points_;
            //this is ownership, not a reference, of course
            mutable autoPtr<edgeList> edgeList_;
            //smallest dp between faces across a fl edge (actually the
            //*largest* dihedral *angle*)
            scalar minDihedral_;
        };

    private:

        const meshRefinement &meshRefiner_;
        const indirectPrimitivePatch &pp_;
        const boolList &regionSnap_;

        //- Internal data
        PtrList<edgeMesh> regionFeatures_;

        DynamicList<featureLine> featureLines_;

        //- Arc length coefficients for all feature lines, populated on demand
        mutable List<List<scalar>* > flParams_;

        //- Search tress for all the feature lines, populated on demand
        mutable List<indexedOctree<treeDataEdge>*> flTrees_;

        //- Search tree for the patch, populated on demand
        mutable autoPtr<indexedOctree<treeDataFace>> patchTree_;

        //- How do we want to extract feature lines?
        bool regionFLs_;
        bool geometryFLs_;

        //- Debug configuration
        bool verbose_; //trigger lots of messages from the chaining algorithm
        bool dump_; //dump vtk files of the features lines and shadow chains

    public:

        featureLinePrep
        (
            const meshRefinement& meshRefiner,
            const indirectPrimitivePatch& pp,
            const dictionary& snapDict,
            const boolList& regionSnap
        );

        ~featureLinePrep();

        //- Main algorithm: constructs feature lines
        void prepareFeatureLines
        (
            const refinementFeatures& rFeatures,
            PtrList<Tuple2<labelList, label>>& featureMeshes
         );

        //- Access to stored fl data
        label getNumFL() const { return featureLines_.size(); }

        //- Access to constructed feature line
        const featureLine &getFeatureLine(label flI) const
        {
            return featureLines_[flI];
        }

        //- Calculation-on-demand of arc length distribution of given
        //feature line
        const List<scalar> &getFLParams(label flI) const;

        //- Calculation-on-demand of edge search tree of given feature line
        const indexedOctree<treeDataEdge> &getFLSearchTree(label flI) const;

    private:

        //- Accessor for the patch tree which builds it if necessary
        const indexedOctree<treeDataFace> &getPatchTree() const;

        //- Converts featureEdgeMesh objects constructed by the feature
        //- refinement algorithm into feature lines for use by the feature
        //- line snapping.
        void mapFeatureMeshes
        (
            const refinementFeatures& rFeatures,
            PtrList<Tuple2<labelList, label>>& featureMeshes
        );

        //- constructs nodes lying between the same two boundaries as sets of
        //- smooth lines
        void extractRegionBdries();

        //- attempts to build one "branch" of a feature line away from the
        //- supplied startEdge in the direction indicated by <reversed> from the
        //- edges in edgeMesh
        void searchFromEdge
        (
            const edgeMesh &edgeMesh,
            PackedList<1> &used,
            label startEdge,
            bool reversed,
            DynamicList<label> &currentFL,
            const labelList &startNodeStartPoints,
            const labelList &endNodeStartPoints,
            const labelList &edgesByStartPoints,
            const labelList &edgesByEndPoints
         ) const;

        //- set up the search tree for a given feature line
        indexedOctree<treeDataEdge>*
        createFeatureLineSearchTree(label flIndex) const;
    };
}
#endif
