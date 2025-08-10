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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include "include/demandDrivenData.H"
#include "utilities/surfaceTools/meshSurfaceEdgeExtractorFUN/meshSurfaceEdgeExtractorFUN.H"
#include "utilities/octrees/meshOctree/meshOctree.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapper.H"
#include "utilities/helperFunctions/helperFunctions.H"
#include "utilities/surfaceTools/createFundamentalSheets/createFundamentalSheetsJFS/createFundamentalSheetsJFS.H"
#include "utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizer.H"
#include "utilities/surfaceTools/meshSurfaceCheckEdgeTypes/meshSurfaceCheckEdgeTypes.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEdgeExtractorFUN::distributeBoundaryFaces()
{
    meshSurfaceEngine mse(mesh_);

    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwner = mse.faceOwners();
    const pointFieldPMG& points = mse.points();

    //- set size of patchNames, newBoundaryFaces_ and newBoundaryOwners_
    const triSurf& surface = meshOctree_.surface();
    const label nPatches = surface.patches().size();

    wordList patchNames(nPatches);
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners(bFaces.size());
    labelLongList newBoundaryPatches(bFaces.size());

    //- set patchNames
    forAll(surface.patches(), patchI)
        patchNames[patchI] = surface.patches()[patchI].name();

    //- append boundary faces
    forAll(bFaces, bfI)
    {
        newBoundaryFaces.appendList(bFaces[bfI]);
        newBoundaryOwners[bfI] = faceOwner[bfI];
    }

    //- find the region for face by finding the patch nearest
    //- to the face centre
    # ifdef USE_OMP
    # pragma omp parallel for if (bFaces.size() > 100) schedule(guided)
    # endif
    forAll(bFaces, bfI)
    {
        const point c = bFaces[bfI].centre(points);

        label facePatch, nt;
        point p;
        scalar distSq;

        meshOctree_.findNearestSurfacePoint(p, distSq, nt, facePatch, c);

        if ((facePatch > -1) && (facePatch < nPatches))
        {
            newBoundaryPatches[bfI] = facePatch;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot distribute a face " << bFaces[bfI] << " into any "
                << "surface patch!. Exiting.." << exit(FatalError);
        }
    }

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );
}

void meshSurfaceEdgeExtractorFUN::reviseCorners()
{

}

void meshSurfaceEdgeExtractorFUN::reviseEdges()
{

}

void meshSurfaceEdgeExtractorFUN::remapBoundaryPoints()
{
    meshSurfaceEngine& mse = surfaceEngine();
    meshSurfaceMapper mapper(mse, meshOctree_);

    mapper.mapVerticesOntoSurfacePatches();
}

void meshSurfaceEdgeExtractorFUN::createBasicFundamentalSheets()
{
    createFundamentalSheetsJFS edgeSheets(mesh_, createWrapperSheet_);

    clearOut();
}

void meshSurfaceEdgeExtractorFUN::smoothMeshSurface()
{
    meshSurfaceEngine& mse = surfaceEngine();

    meshSurfaceOptimizer mso(mse, meshOctree_);
    mso.optimizeSurface();
}

void meshSurfaceEdgeExtractorFUN::improveQualityOfFundamentalSheets()
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const edgeList& edges = mse.edges();

    meshSurfaceCheckEdgeTypes edgeCheck(mse);

    label id = mesh_.addPointSubset("convexEdges");
    labelLongList helper;
    edgeCheck.convexEdges(helper);
    forAll(helper, i)
    {
        mesh_.addPointToSubset(id, edges[helper[i]].start());
        mesh_.addPointToSubset(id, edges[helper[i]].end());
    }

    id = mesh_.addPointSubset("concaveEdges");
    edgeCheck.concaveEdges(helper);
    forAll(helper, i)
    {
        mesh_.addPointToSubset(id, edges[helper[i]].start());
        mesh_.addPointToSubset(id, edges[helper[i]].end());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
