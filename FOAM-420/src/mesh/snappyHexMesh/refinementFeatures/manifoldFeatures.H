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
    (c) 2017 Esi Ltd.

Class
    Foam::manifoldFeatures

Description
    Construct manifold feature edges together with seeds containing
    mesh cell label closest to feature start. Used for seeding and accelerating
    feature refinement tracking

SourceFiles
    manifoldFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef manifoldFeatures_H
#define manifoldFeatures_H

#include "edgeMesh/extendedEdgeMesh/extendedFeatureEdgeMesh/extendedFeatureEdgeMesh.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "primitives/Tuple2/Tuple2.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class manifoldFeatures Declaration
\*---------------------------------------------------------------------------*/

class manifoldFeatures
{
private:

    //- List of extended edge meshes
    const PtrList<extendedFeatureEdgeMesh>& extFeatures_;

    //- Per shell if refinement feature only
    const List<bool>& refinementOnly_;

    //- Whether to check for refinement only features
    const bool& checkRefinementOnly_;

    //- Manifold feature edges
    PtrList<Tuple2<labelList, label>> features_;

    //- Closest mesh cell to starting point
    labelList seeds_;

public:

    //- Runtime type information
    ClassName("manifoldFeatures");

    // Constructors

    manifoldFeatures
    (
        const PtrList<extendedFeatureEdgeMesh>& extFeatures,
        const List<bool>& refinementOnly,
        const bool& checkRefinementOnly
    )
    :
        extFeatures_(extFeatures),
        refinementOnly_(refinementOnly),
        checkRefinementOnly_(checkRefinementOnly),
        features_(0),
        seeds_(0)
    {}

    //- Return a list of manifold features
    PtrList<Tuple2<labelList, label>>& features
    (
        const label maxFtrSz = labelMax,
        const scalar minFtrLength = scalar(-1)
    );

    //- Constant access to seed cells
    const labelList& seeds() const
    {
        return seeds_;
    }

    //- Non-constant access to seed cells
    labelList& seeds()
    {
        return seeds_;
    }

    //- Clear manifold edges and seed cells
    void clear();

    //- Recalculate seed cells on re-distrubution
    void distribute
    (
        const mapDistributePolyMesh& map
    );

    //- Recalculate seed cells on mesh update
    void update
    (
        const mapPolyMesh& map
    );

    //- Calculate feature length for filtering
    scalar calcFeatureLength
    (
        const edgeMesh& eMesh,
        const labelList& meshPoints
    );

    //- Calculate a set of manifold edges
    void calcManifoldFeatureMeshes
    (
        PtrList<Tuple2<labelList,label>>& manifoldFeatures,
        const label maxFtrSz = labelMax,
        const scalar minFtrLength = scalar(-1)
    );

    //- Perform walk along feature edge mesh from seed point
    static inline void manifoldEdgeWave
    (
        const label& prevPointI,
        const label& pointI,
        const labelListList& pointEdges,
        const edgeMesh& featureMesh,
        PackedList<1>& edgeVisited,
        PackedList<1>& pointVisited,
        DynamicList<label>& meshPoints,
        const label maxFtrSz
    )
    {
        SortableList<scalar> maxDot(pointEdges[pointI].size());
        const pointField& points = featureMesh.points();
        vector vec = points[pointI] - points[prevPointI];
        vec /= (mag(vec) + SMALL);

        forAll(pointEdges[pointI], i)
        {
            label edgeI = pointEdges[pointI][i];
            const edge& e = featureMesh.edges()[edgeI];
            label nextPointI = e.otherVertex(pointI);
            if (nextPointI != prevPointI)
            {
                vector vec1 = points[nextPointI] - points[pointI];
                vec1 /= (mag(vec1) + SMALL);
                maxDot[i] = -(vec & vec1);
            }
            else
            {
                maxDot[i] = 1;
            }
        }
        maxDot.sort();

        const labelList& indices = maxDot.indices();

        forAll(indices, i)
        {
            label index = indices[i];
            label edgeI = pointEdges[pointI][index];
            const edge& e = featureMesh.edges()[edgeI];
            label nextPointI = e.otherVertex(pointI);

            if
            (
                (
                    edgeVisited.get(edgeI) == 0 && maxDot[i] > -0.01//0.99
                    && maxDot.size() > 1
                )
                || meshPoints.size() > maxFtrSz
            )
            {
                return;
            }

            if
            (
                edgeVisited.set(edgeI, 1u)
            )
            {
                meshPoints.append(nextPointI);

                if
                (
                    pointVisited.set(nextPointI, 1u)
                )
                {
                    manifoldEdgeWave
                    (
                        pointI,
                        nextPointI,
                        pointEdges,
                        featureMesh,
                        edgeVisited,
                        pointVisited,
                        meshPoints,
                        maxFtrSz
                    );
                }
                break;
            }
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
