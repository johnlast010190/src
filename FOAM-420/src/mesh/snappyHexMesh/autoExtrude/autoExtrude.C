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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "autoExtrude/autoExtrude.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "shellSurfaces/shellSurfaces.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "polyTopoChange/polyTopoChange/addPatchCellLayer.H"
#include "autoOptimize/autoOptimize.H"
#include "edgeClassification/edgeClassification.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoExtrude, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::autoExtrude::extType,
        4
    >::names[] =
    {
        "direction",
        "patchNormal",
        "avePatchNormal",
        "target"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::autoExtrude::extAction,
        3
    >::names[] =
    {
        "add",
        "new",
        "split"
    };
}

const Foam::NamedEnum<Foam::autoExtrude::extType, 4>
    Foam::autoExtrude::extTypeNames;

const Foam::NamedEnum<Foam::autoExtrude::extAction, 3>
    Foam::autoExtrude::extActionNames;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::autoExtrude::autoExtrude
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::autoExtrude::preBalance
(
    const dictionary& extrudeDict,
    const bool updateSurf
)
{
    Info<<"Pre-balance before extruding"<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    scalarField cellWeights(mesh.nCells(), scalar(1));
    forAllConstIter(dictionary, extrudeDict, iter)
    {
        const word& key = iter().keyword();
        const dictionary& edict = extrudeDict.subDict(key);

        if (!edict.found("sourcePatches"))
        {
            continue;
        }

        label nLayers = edict.lookupOrDefault<label>("nLayers", 1);

        labelHashSet sourceSet = pbm.patchSet
        (
            wordReList(edict.lookup("sourcePatches")),
            false,
            true
        );

        filterExtrudePatch(sourceSet.toc());
        autoPtr<indirectPrimitivePatch> extrudePatchPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                sourceSet.toc()
            )
        );
        const indirectPrimitivePatch& extrudePatch = extrudePatchPtr();

        forAll(extrudePatch, i)
        {
            label facei = extrudePatch.addressing()[i];
            label own = mesh.faceOwner()[facei];
            cellWeights[own] += scalar(nLayers);
        }
    }

    meshRefiner_.balance
    (
        false,
        true,    // keepZoneFaces
        false,   // keepBaffles
        cellWeights,
        decomposer_,
        distributor_,
        updateSurf
    );

    return;
}


void Foam::autoExtrude::balance
(
    const scalar snapWeights,
    const bool updateSurf
)
{
    Info<<"Rebalance extruded mesh"<<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    scalarField weightField(mesh.nCells(),1);
    if (snapWeights > 0)
    {
        meshRefiner_.calcSnapWeights
        (
            snapWeights,
            weightField
        );
    }

    meshRefiner_.balance
    (
        false,
        true,    // keepZoneFaces
        false,   // keepBaffles
        weightField,//scalarField(mesh.nCells(), 1),
        decomposer_,
        distributor_,
        updateSurf
    );

    return;
}


void Foam::autoExtrude::optimize
(
    const refinementParameters& refineParams
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //Inflate by smoothing internal points
    labelList boundaryPts(mesh.nPoints(), -1);
    forAll(mesh.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);
        if (patchi != -1 && !patches[patchi].coupled())
        {
            const face& f = mesh.faces()[facei];
            forAll(f, fp)
            {
                boundaryPts[f[fp]] = patchi;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPts,
        maxEqOp<label>(),
        label(-1)               // null value
    );

    DynamicList<label> bPatches(patches.size());
    forAll(patches ,patchi)
    {
        if (!patches[patchi].coupled())
        {
            bPatches.append(patchi);
        }
    }

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            bPatches
         )
    );
    indirectPrimitivePatch& pp = ppPtr();
    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    boolList excludedFaces(pp.size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        pp,
        meshEdges,
        excludedFaces,
        0.707,
        0.707
    );

    pointField surfNormals = eClass.calculatePointNormals
    (
        excludedFaces,
        0,
        true
    );

    boolList layerEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        const edge& e = mesh.edges()[edgei];
        if
        (
            (boundaryPts[e[0]] != -1 && boundaryPts[e[1]] == -1)
            || (boundaryPts[e[1]] != -1 && boundaryPts[e[0]] == -1)
        )
        {
            layerEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    labelList maxPointLevel(mesh.nPoints(), -labelMax);

    forAll(cellLevel, celli)
    {
        const labelList& cPts = mesh.cellPoints()[celli];
        label cLevel = cellLevel[celli];

        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];

            maxPointLevel[pointi] = max(maxPointLevel[pointi],cLevel);
        }
    }
    syncTools::syncPointList
    (
        mesh,
        maxPointLevel,
        maxEqOp<label>(),
        label(-1)
     );

    vectorField dispVec(mesh.nPoints(), vector::zero);
    forAll(surfNormals, pti)
    {
        label meshpointi = meshPoints[pti];
        const labelList& pEdges = mesh.pointEdges()[meshpointi];
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (layerEdges[edgei])
            {
                edge e = mesh.edges()[edgei];
                label otherPt
                (
                    e[0] == meshpointi ? e[1] : e[0]
                );

                vector sN = -surfNormals[pti];

                point startPt = mesh.points()[meshpointi];
                point outerPt = mesh.points()[otherPt];
                scalar edgeLen = edge0Len /
                    (1<<maxPointLevel[meshpointi]);

                point projPt = startPt + 0.33*edgeLen*sN;
                dispVec[otherPt] =  projPt - outerPt;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        dispVec,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    dispVec += mesh.points();
    mesh.movePoints(dispVec);

    // Optimise and inflate layers using cell sphericity
    dictionary dummyOptimDict;
    dummyOptimDict.add("type","foamOptimize",true);
    dummyOptimDict.add
        ("foamOptimizeCoeffs",dictionary(), true);
    dummyOptimDict.subDict("foamOptimizeCoeffs").add
    (
        "optimizationMethod",
        "fullMeshOptimization",
        true
    );
    dummyOptimDict.subDict("foamOptimizeCoeffs").add
    (
        "anisotropy",
        "isotropy",
        true
    );
    dummyOptimDict.subDict("foamOptimizeCoeffs").add
    (
        "iterations",
        "6",
        true
    );

    autoPtr<autoOptimize> optimMeshPtr
        = autoOptimize::New(mesh, dummyOptimDict);
    optimMeshPtr->optimize();

    //Move cells on wrong side of geometry
    meshRefiner_.moveWrongSidedCells(refineParams,false);
    optimMeshPtr->optimize();

    return;
}


bool Foam::autoExtrude::calculateDirectionDisplacement
(
    const dictionary& edict,
    const indirectPrimitivePatch& extrudePatch,
    const regionSplit2D& sourceRegionIDs,
    pointField& displacement
)
{
    scalarField sourceExtrudeFCH;
    pointField sourceRegionCentroids;
    pointField sourceRegionNormals;
    scalarField sourceRegionAreas;
    calcRegionsProperties
    (
        extrudePatch,
        sourceRegionIDs,
        sourceExtrudeFCH,
        sourceRegionCentroids,
        sourceRegionNormals,
        sourceRegionAreas
    );

    vector dir = edict.lookup("direction");
    vector uDir = dir/mag(dir);
    if (edict.found("fch"))
    {
        scalar extrudeFCH = readScalar(edict.lookup("fch"));
        displacement = uDir*extrudeFCH;
    }
    else
    {
        forAll(extrudePatch, i)
        {
            label regioni = sourceRegionIDs[i];
            scalar extrudeFCH = sourceExtrudeFCH[regioni];
            const face& f =  extrudePatch.localFaces()[i];
            forAll(f,fp)
            {
                displacement[f[fp]] = uDir*extrudeFCH;
            }
        }
    }
    return true;
}


bool Foam::autoExtrude::calculateNormalDisplacement
(
    const dictionary& edict,
    const indirectPrimitivePatch& extrudePatch,
    const regionSplit2D& sourceRegionIDs,
    pointField& displacement
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    scalarField sourceExtrudeFCH;
    pointField sourceRegionCentroids;
    pointField sourceRegionNormals;
    scalarField sourceRegionAreas;
    calcRegionsProperties
    (
        extrudePatch,
        sourceRegionIDs,
        sourceExtrudeFCH,
        sourceRegionCentroids,
        sourceRegionNormals,
        sourceRegionAreas
    );

    displacement = PatchTools::pointNormals(mesh, extrudePatch);
    vectorField unitDisp(displacement.size(), vector::zero);
    unitDisp = displacement/mag(displacement);

    if (edict.found("fch"))
    {
        scalar extrudeFCH  = readScalar(edict.lookup("fch"));
        displacement = unitDisp*extrudeFCH;
    }
    else
    {
        forAll(extrudePatch, i)
        {
            label regioni = sourceRegionIDs[i];
            scalar extrudeFCH = sourceExtrudeFCH[regioni];
            const face& f =  extrudePatch.localFaces()[i];
            forAll(f,fp)
            {
                displacement[f[fp]] = unitDisp[f[fp]]*extrudeFCH;
            }
        }
    }

    return true;
}


bool Foam::autoExtrude::calculateAverageNormalDisplacement
(
    const dictionary& edict,
    const indirectPrimitivePatch& extrudePatch,
    const regionSplit2D& sourceRegionIDs,
    pointField& displacement
)
{
    scalarField sourceExtrudeFCH;
    pointField sourceRegionCentroids;
    pointField sourceRegionNormals;
    scalarField sourceRegionAreas;
    calcRegionsProperties
    (
        extrudePatch,
        sourceRegionIDs,
        sourceExtrudeFCH,
        sourceRegionCentroids,
        sourceRegionNormals,
        sourceRegionAreas
    );

    scalar extrudeFCH = -1;
    bool useInputFCH = false;
    if (edict.found("fch"))
    {
        useInputFCH = true;
        extrudeFCH  = readScalar(edict.lookup("fch"));
    }

    forAll(extrudePatch, i)
    {
        label regioni = sourceRegionIDs[i];
        if (!useInputFCH)
        {
            extrudeFCH = sourceExtrudeFCH[regioni];
        }

        vector extrudeDir = sourceRegionNormals[regioni];
        vector unitDir = extrudeDir / mag(extrudeDir);
        const face& f =  extrudePatch.localFaces()[i];
        forAll(f,fp)
        {
            displacement[f[fp]] = extrudeFCH*unitDir;
        }
    }

    return true;
}


Foam::scalar Foam::autoExtrude::calculateStretch
(
    const scalar lheight,
    const scalar fch,
    const label nLayers
)
{
    bool oddNLayers((nLayers % 2) == 0 ? false : true);
    label halfLayers = nLayers/2;
    scalar halfHeight = lheight/2;
    scalar tol = 1e-6;
    scalar xn = 1.1;

    //Guess a starting stretch
    while (xn < 5)
    {
        scalar h = fch*(1.- pow(xn,halfLayers))/(1.-xn);
        if (h > halfHeight)
        {
            break;
        }
        xn *= 1.1;
    }

    scalar x = -1.0;
    while (mag(x - xn) > tol)
    {
        x = xn;
        scalar f1 = lheight-2.0*fch;
        scalar f2 = 0;
        for (int i = 1; i < halfLayers; i++)
        {
            f1 -= 2.0*fch*pow(x,i);
            f2 -= 2.0*i*fch*pow(x,i-1);
        }
        if (oddNLayers)
        {
            f1 -= fch*pow(x,halfLayers);
            f2 -= (halfLayers)*fch*pow(x,halfLayers-1);
        }
        xn = x - f1/(f2+SMALL);
    }

    if (xn < 1.0 + tol && xn > 1.0 - tol)
    {
        xn = 0.0;
        x = -1.0;
        while (mag(x - xn) > tol)
        {
            x = xn;
            scalar f1 = lheight-2.0*fch;
            scalar f2 = 0;
            for (int i = 1; i < halfLayers; i++)
            {
                f1 -= 2.0*fch*pow(x,i);
                f2 -= 2.0*i*fch*pow(x,i-1);
            }
            if (oddNLayers)
            {
                f1 -= fch*pow(x,halfLayers);
                f2 -= (halfLayers)*fch*pow(x,halfLayers-1);
            }
            xn = x - f1/(f2+SMALL);
        }
    }

    return xn;
}

bool Foam::autoExtrude::calculateTargetDisplacement
(
    const dictionary& edict,
    const indirectPrimitivePatch& extrudePatch,
    const labelHashSet& sourceSet,
    const regionSplit2D& sourceRegionIDs,
    const label nLayers,
    pointField& displacement,
    scalarField& expansionRatio
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    scalarField sourceExtrudeFCH;
    pointField sourceRegionCentroids;
    pointField sourceRegionNormals;
    scalarField sourceRegionAreas;
    calcRegionsProperties
    (
        extrudePatch,
        sourceRegionIDs,
        sourceExtrudeFCH,
        sourceRegionCentroids,
        sourceRegionNormals,
        sourceRegionAreas
     );

    Random rndGen(653213);

    labelHashSet targetSet = pbm.patchSet
    (
        wordReList(edict.lookup("targetPatches")),
        false,
        true
    );

    autoPtr<indirectPrimitivePatch> targetPatchPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            targetSet.toc()
         )
    );
    const indirectPrimitivePatch& targetPatch = targetPatchPtr();

    Switch flatten = edict.lookupOrDefault<Switch>("flatten", false);
    if (flatten)
    {
        flattenPatch(targetPatch);
    }

    if (returnReduce(targetPatch.size(), sumOp<label>()) == 0)
    {
        WarningInFunction
            <<"Extrusion contains zero sized target patch,"
            <<" skipping entry : " << edict.name()
            << nl << endl;
        return false;
    }

    boolList blockedTargetFaces(targetPatch.size(), true);
    regionSplit2D targetRegionIDs(mesh, targetPatch, blockedTargetFaces);

    label nTargetRegions = targetRegionIDs.nRegions();
    label nSourceRegions = sourceRegionIDs.nRegions();

    if (nSourceRegions == nTargetRegions)
    {
        Switch translational =
            edict.lookupOrDefault<Switch>("translational", true);
        if (!translational)
        {
            // Precalculate mesh edges for pp.edges.
            const labelList meshEdges
            (
                extrudePatch.meshEdges
                (
                    mesh.edges(),
                    mesh.pointEdges()
                 )
             );
            labelList nExternalEdge(mesh.nEdges(), 0);
            labelList sourceRegionEdgeID(mesh.nEdges(), -1);

            const labelListList& edgeFaces = extrudePatch.edgeFaces();

            forAll(extrudePatch.edges(), edgeI)
            {
                label meshEdgeI = meshEdges[edgeI];
                const labelList& eFaces = edgeFaces[edgeI];
                forAll(eFaces, eFI)
                {
                    sourceRegionEdgeID[meshEdgeI] = max
                    (
                        sourceRegionIDs[eFaces[eFI]],
                        sourceRegionEdgeID[meshEdgeI]
                    );
                }
                nExternalEdge[meshEdgeI] = eFaces.size();
            }

            syncTools::syncEdgeList
            (
                mesh,
                nExternalEdge,
                plusEqOp<label>(),
                label(0)              // null value
            );
            syncTools::syncEdgeList
            (
                mesh,
                sourceRegionEdgeID,
                maxEqOp<label>(),
                label(-1)              // null value
             );

            vectorField sideDir(nSourceRegions, vector::zero);
            labelList nSideEdges(nSourceRegions, label(0));

            forAll(mesh.edges(), edgei)
            {
                if (nExternalEdge[edgei] == 1)
                {
                    label regioni = sourceRegionEdgeID[edgei];
                    const labelList& eFaces =  mesh.edgeFaces()[edgei];
                    forAll(eFaces, eFI)
                    {
                        label facei = eFaces[eFI];
                        label patchi = pbm.whichPatch(facei);
                        if (patchi != -1 && !pbm[patchi].coupled())
                        {
                            if (!sourceSet.found(patchi))
                            {
                                const edge& e = mesh.edges()[edgei];
                                vector eVec =  e.vec(mesh.points());
                                vector fVec = mesh.faceAreas()[facei];

                                vector gDir(eVec ^ fVec);
                                scalar magGDir = mag(gDir);
                                if (magGDir > SMALL)
                                {
                                    gDir /= magGDir;
                                    if ((gDir&sourceRegionNormals[regioni]) < 0)
                                    {
                                        gDir = -gDir;
                                    }
                                    sideDir[regioni] += gDir;
                                    nSideEdges[regioni]++;
                                }
                            }
                        }
                    }
                }
            }

            Pstream::listCombineReduce(sideDir, plusOp<vector>());
            Pstream::listCombineReduce(nSideEdges, plusOp<label>());

            forAll(sideDir, regioni)
            {
                if (nSideEdges[regioni] > 0)
                {
                    sideDir[regioni] /= nSideEdges[regioni];
                    sideDir[regioni] /= mag(sideDir[regioni]);
                }
            }
            sourceRegionNormals = sideDir;
        }

        const pointField& mfCtrs = extrudePatch.faceCentres();
        pointField targetRegionCentroids(nTargetRegions,vector::zero);
        pointField targetRegionNormals(nTargetRegions,vector::zero);
        scalarField targetRegionAreas(nTargetRegions,scalar(0));

        scalarField sourceRegionRadius(nSourceRegions,scalar(0));

        forAll(sourceRegionIDs, facei)
        {
            label regioni = sourceRegionIDs[facei];
            scalar radius = mag
            (
                mfCtrs[facei] - sourceRegionCentroids[regioni]
            );
            sourceRegionRadius[regioni] = max
            (
                radius,
                sourceRegionRadius[regioni]
            );
        }

        Pstream::listCombineReduce(sourceRegionRadius, maxOp<scalar>());

        const pointField& sfCtrs = targetPatch.faceCentres();
        const vectorField& sfAreas = targetPatch.faceAreas();
        const vectorField& sfNormals = targetPatch.faceNormals();

        forAll(targetRegionIDs, facei)
        {
            label regioni = targetRegionIDs[facei];
            scalar fArea = mag(sfAreas[facei]);
            targetRegionCentroids[regioni] += sfCtrs[facei]*fArea;
            targetRegionNormals[regioni] += sfNormals[facei]*fArea;
            targetRegionAreas[regioni] += fArea;
        }

        Pstream::listCombineReduce(targetRegionCentroids, plusOp<vector>());
        Pstream::listCombineReduce(targetRegionNormals, plusOp<vector>());
        Pstream::listCombineReduce(targetRegionAreas, plusOp<scalar>());

        forAll(targetRegionCentroids, regioni)
        {
            targetRegionCentroids[regioni]
                /= targetRegionAreas[regioni];
            targetRegionNormals[regioni]
                /= targetRegionAreas[regioni];
        }

        boundBox sourceTargetBB(extrudePatch.localPoints());
        boundBox targetBB(targetPatch.localPoints());
        sourceTargetBB.add(targetBB);
        treeBoundBox treeCombinedBb(sourceTargetBB);
        treeCombinedBb = treeCombinedBb.extend(rndGen, 1E-4);

        indexedOctree<treeDataFace> targetTree
        (
            treeDataFace(false, mesh, targetPatch.addressing()),
            treeCombinedBb,  // overall search domain
            12,         // maxLevel
            10,         // leafsize
            6.0         // duplicity
         );

        scalar bbSize = treeCombinedBb.mag();

        labelList targetHitRegions(nSourceRegions, -1);
        for (label regioni=0; regioni < nSourceRegions; ++regioni)
        {
            scalar tol = 0.01*sourceRegionRadius[regioni];
            for (int i = 0; i < 2; i++)
            {
                for (direction j = 0; j < vector::nComponents; j++)
                {
                    vector dir = vector::zero;
                    if (i == 0)
                    {
                        dir[j] = tol;
                    }
                    else
                    {
                        dir[j] = -tol;
                    }

                    if (targetHitRegions[regioni] == -1)
                    {
                        point start =
                            sourceRegionCentroids[regioni] + dir;
                        vector mNormal = sourceRegionNormals[regioni];
                        mNormal /= mag(mNormal);
                        point end  = start + bbSize*mNormal;
                        pointIndexHit hitInfo =
                            targetTree.findLine(start,end);
                        if (hitInfo.hit())
                        {
                            label index = hitInfo.index();
                            targetHitRegions[regioni] =
                                targetRegionIDs[index];
                        }
                    }
                }
            }
        }
        Pstream::listCombineReduce(targetHitRegions, maxOp<label>());

        //If any targets unset search for a nearest one
        bool unset = false;
        for (label regioni=0; regioni < nSourceRegions; ++regioni)
        {
            if (targetHitRegions[regioni] == -1)
            {
                point sCentroid = sourceRegionCentroids[regioni];
                vector mNormal = sourceRegionNormals[regioni];
                mNormal /= mag(mNormal);
                scalar maxOrient = -GREAT;
                label tRegion = -1;
                for (label tregioni=0; tregioni < nTargetRegions; ++tregioni)
                {
                    point tCentroid = targetRegionCentroids[tregioni];
                    vector sourceToTarget = (tCentroid-sCentroid);
                    scalar sourceToTargetDistance = mag(sourceToTarget);
                    if (sourceToTargetDistance > 0)
                    {
                        sourceToTarget /= sourceToTargetDistance;
                        scalar orient(mag(sourceToTarget&mNormal));
                        if (orient >  maxOrient)
                        {
                            tRegion = tregioni;
                            maxOrient = orient;
                        }
                    }
                }
                WarningInFunction
                    << "Unable to find target for source region " << regioni
                    << " with centroid : " << sCentroid
                    << " Matching with centroid : "
                    << targetRegionCentroids[tRegion]<<endl;
                targetHitRegions[regioni] = tRegion;
                unset = true;
            }
        }
        if (returnReduce(unset, orOp<bool>()))
        {
            Pstream::listCombineReduce(targetHitRegions, maxOp<label>());
        }

        if (translational)
        {
            scalar matchingAreaTol = 0.95;
            vectorField sourceExtrudeDir(nSourceRegions, vector::zero);
            forAll(sourceExtrudeDir, regioni)
            {
                label targetRegion = targetHitRegions[regioni];
                point targetCentre = targetRegionCentroids[targetRegion];
                point sourceCentre = sourceRegionCentroids[regioni];

                scalar sArea = sourceRegionAreas[targetRegion];
                scalar tArea = targetRegionAreas[targetRegion];

                scalar aRatio = min(sArea,tArea)/max(sArea,tArea);

                if (aRatio < matchingAreaTol)
                {
                    WarningInFunction
                        <<"Matching area ratio : " << aRatio
                        <<" less than tolerance : " << matchingAreaTol
                        <<" for source region centre " << sourceCentre
                        <<" for target region centre " << targetCentre
                        <<" continuing but this might be a problem."
                        << nl << endl;
                }

                sourceExtrudeDir[regioni] = targetCentre - sourceCentre;
            }

            scalarField sourceExtrudeExpansion(nSourceRegions, scalar(1.0));
            forAll(sourceExtrudeExpansion, regioni)
            {
                scalar lheight = mag(sourceExtrudeDir[regioni]);
                scalar fch = sourceExtrudeFCH[regioni];
                //Calculate expansion ratio
                sourceExtrudeExpansion[regioni] = calculateStretch
                (
                    lheight,
                    fch,
                    nLayers
                );
            }

            //Set displacement and stretching
            forAll(extrudePatch, i)
            {
                label regioni = sourceRegionIDs[i];
                scalar extrudeFCH = sourceExtrudeFCH[regioni];
                vector extrudeDir = sourceExtrudeDir[regioni];
                vector unitDir = extrudeDir / mag(extrudeDir);

                const face& f =  extrudePatch.localFaces()[i];
                forAll(f,fp)
                {
                    displacement[f[fp]] = extrudeFCH*unitDir;
                    expansionRatio[f[fp]] = sourceExtrudeExpansion[regioni];
                }
            }
        }
        else
        {
            const labelList& sourceMeshPts = extrudePatch.meshPoints();
            labelList sourcePtRegionID(sourceMeshPts.size(), -1);
            forAll(extrudePatch, i)
            {
                label regioni = sourceRegionIDs[i];
                const face& f =  extrudePatch.localFaces()[i];
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    sourcePtRegionID[pointi] = max
                    (
                        regioni,
                        sourcePtRegionID[pointi]
                     );
                }
            }

            forAll(sourceMeshPts, ptI)
            {
                label regioni = sourcePtRegionID[ptI];
                label meshpointi = sourceMeshPts[ptI];

                point sourceCentre = sourceRegionCentroids[regioni];

                label targetRegion = targetHitRegions[regioni];
                point targetCentre = targetRegionCentroids[targetRegion];
                point targetNormal = targetRegionNormals[targetRegion];
                targetNormal /= mag(targetNormal);

                vector targetDir = targetCentre - sourceCentre;

                point startPt = mesh.points()[meshpointi];
                point endPt = startPt + 2*targetDir;

                linePointRef targetLine(startPt, endPt);
                plane tpl(targetCentre, targetNormal, false);

                scalar cutPt = tpl.lineIntersect(targetLine);
                scalar layerHeight = mag(cutPt*targetLine.vec());

                scalar fch = sourceExtrudeFCH[regioni];
                scalar str = calculateStretch
                (
                    layerHeight,
                    fch,
                    nLayers
                );

                vector unitDir = targetDir / mag(targetDir);
                displacement[ptI] = fch*unitDir;
                expansionRatio[ptI] = str;
            }
        }
    }
    else
    {
        WarningInFunction
            <<"Number of source regions : " << nSourceRegions
            <<" does not match number of target regions  : "
            << nTargetRegions <<". Skipping entry : " << edict.name()
            << nl << endl;
        return false;
    }

    return true;
}


void Foam::autoExtrude::removeUnusedPoints()
{
    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);

    label nPointsRemoved = 0;
    boolList removedPts(mesh.nPoints(), false);

    const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
    labelList nPointEdges(mesh.nPoints(), 0);

    forAll(mesh.points(), pointI)
    {
        const labelList pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            label edgeI = pEdges[pEI];
            if (isMasterEdge[edgeI])
            {
                nPointEdges[pointI]++;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        nPointEdges,
        plusEqOp<label>(),
        label(0)
    );

    labelList nPointFaces(mesh.nPoints(), 0);
    forAll(mesh.points(), pointI)
    {
        const labelList pFaces = mesh.pointFaces()[pointI];
        forAll(pFaces, pFI)
        {
            label faceI = pFaces[pFI];
            if (isMasterFace[faceI])
            {
                nPointFaces[pointI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nPointFaces,
        plusEqOp<label>(),
        label(0)
    );


    forAll(mesh.faces(), facei)
    {
        face f = mesh.faces()[facei];

        DynamicList<label> newFacePts(f.size());

        forAll(f, fp)
        {
            if (nPointEdges[f[fp]] == 2 && nPointFaces[f[fp]] > 2)
            {
                if (!removedPts[f[fp]])
                {
                    meshMod.setAction(polyRemovePoint(f[fp]));
                    removedPts[f[fp]] = true;
                    nPointsRemoved++;
                }
            }
            else
            {
                newFacePts.append(f[fp]);
            }
        }
        newFacePts.shrink();

        face newFace(newFacePts);
        if (newFace.size() != f.size())
        {
            label patchID = mesh.boundaryMesh().whichPatch(facei);
            label nei =
                (patchID == -1 ? mesh.faceNeighbour()[facei] : -1);

            label zoneID = mesh.faceZones().whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    face(newFace),              // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],    // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchID,                    // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                 )
             );
        }
    }

    if (returnReduce(nPointsRemoved, sumOp<label>()))
    {
        // Create mesh (no inflation), return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Optionally inflate mesh
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh.setInstance(meshRefiner_.timeName());

        // Update intersection info
        meshRefiner_.updateMesh(map, labelList(0));
    }
}


void Foam::autoExtrude::addNewProcessorPatches
(
    const indirectPrimitivePatch& pp,
    labelList& edgePatchID,
    labelList& edgeZoneID,
    boolList& edgeFlip,
    labelList& inflateFaceID
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
         )
    );

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches)
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
         )
    );

    label nPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;
    addPatchCellLayer::calcExtrudeInfo
    (
        true,   // for internal edges get zone info from any face

        mesh,
        globalFaces,
        edgeGlobalFaces,
        pp,
        labelList(0), //list of grown-up patches

        edgePatchID,
        nPatches,
        nbrProcToPatch,
        patchToNbrProc,

        edgeZoneID,
        edgeFlip,
        inflateFaceID
     );

    label nAdded = nPatches - mesh.boundaryMesh().size();
    reduce(nAdded, sumOp<label>());

    Info<< "Adding overall " << nAdded << " processor patches." << endl;

    if (nAdded > 0)
    {
        DynamicList<polyPatch*> newPatches(nPatches);
        forAll(mesh.boundaryMesh(), patchi)
        {
            newPatches.append
            (
                mesh.boundaryMesh()[patchi].clone
                (
                    mesh.boundaryMesh()
                 ).ptr()
             );
        }

        for
        (
            label patchi = mesh.boundaryMesh().size();
            patchi < nPatches;
            patchi++
        )
        {
            label nbrProci = patchToNbrProc[patchi];
            word name
            (
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProci)
            );

            dictionary patchDict;
            patchDict.add("type", processorPolyPatch::typeName);
            patchDict.add("myProcNo", Pstream::myProcNo());
            patchDict.add("neighbProcNo", nbrProci);
            patchDict.add("nFaces", 0);
            patchDict.add("startFace", mesh.nFaces());

            meshRefiner_.appendPatch
            (
                mesh,
                name,
                patchDict
            );
        }
        mesh.clearOut();
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
    }
}


void Foam::autoExtrude::calcRegionsProperties
(
    const indirectPrimitivePatch& pp,
    const regionSplit2D& regionIDs,

    scalarField& extrudeFCH,
    pointField& regionCentroids,
    pointField& regionNormals,
    scalarField& regionAreas
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    label nRegions = regionIDs.nRegions();
    extrudeFCH.setSize(nRegions, scalar(0));
    regionCentroids.setSize(nRegions,vector::zero);
    regionNormals.setSize(nRegions,vector::zero);
    regionAreas.setSize(nRegions,scalar(0));

    const vectorField& fNormals = pp.faceNormals();
    const pointField& fCtrs = pp.faceCentres();
    const vectorField& fAreas = pp.faceAreas();

    forAll(regionIDs, i)
    {
        label facei = pp.addressing()[i];
        label own = mesh.faceOwner()[facei];
        scalar edgeLen = edge0Len / (1<<cellLevel[own]);
        label regioni = regionIDs[i];
        scalar fArea = mag(fAreas[i]);
        extrudeFCH[regioni] += edgeLen*fArea;
        regionCentroids[regioni] += fCtrs[i]*fArea;
        regionNormals[regioni] += fNormals[i]*fArea;
        regionAreas[regioni] += fArea;
    }

    Pstream::listCombineReduce(extrudeFCH, plusOp<scalar>());
    Pstream::listCombineReduce(regionCentroids, plusOp<vector>());
    Pstream::listCombineReduce(regionNormals, plusOp<vector>());
    Pstream::listCombineReduce(regionAreas, plusOp<scalar>());

    forAll(regionCentroids, regioni)
    {
        extrudeFCH[regioni]
            /= regionAreas[regioni];
        regionCentroids[regioni]
            /= regionAreas[regioni];
        regionNormals[regioni]
            /= regionAreas[regioni];
    }

    return;
}


void Foam::autoExtrude::filterExtrudePatch
(
    const labelList& patchIDs
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList addressing;
    meshRefinement::calcPatchAddressing(mesh,patchIDs,addressing);

    labelHashSet extrudedSet(patchIDs);

    //Block faces leading to bad extrusion
    boolList blockedFaces(mesh.nFaces(), false);
    labelList nBoundaryFaces(mesh.nPoints(), 0);
    labelList nbrPatchID(mesh.nPoints(), -1);
    forAll(mesh.faces(), facei)
    {
        label patchi = pbm.whichPatch(facei);
        if (patchi != -1 && !pbm[patchi].coupled())
        {
            const face& f = mesh.faces()[facei];
            bool foundNbrPatch = false;
            if (!extrudedSet.found(patchi))
            {
                foundNbrPatch = true;
            }

            forAll(f,fp)
            {
                label pointi = f[fp];
                nBoundaryFaces[pointi]++;
                if (foundNbrPatch)
                {
                    nbrPatchID[pointi] = max(patchi, nbrPatchID[pointi]);
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nBoundaryFaces,
        plusEqOp<label>(), // combine op
        label(0)     // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nbrPatchID,
        maxEqOp<label>(), // combine op
        label(-1)     // null value
    );

    polyTopoChange meshMod(mesh);
    label nUpdated = 0;
    forAll(addressing, i)
    {
        label facei = addressing[i];
        const face& f = mesh.faces()[facei];
        if (f.size() == 3)
        {
            label newPatchID = -1;
            bool nonMan = false;
            forAll(f,fp)
            {
                label pointi = f[fp];

                if (nBoundaryFaces[pointi] == 2)
                {
                    nonMan = true;
                }
                newPatchID = max(nbrPatchID[pointi], newPatchID);
            }

            if (nonMan && newPatchID != -1)
            {
                label own = mesh.faceOwner()[facei];
                label nei = -1;

                label zoneI = mesh.faceZones().whichZone(facei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    label index = mesh.faceZones()[zoneI].whichFace(facei);
                    zoneFlip = mesh.faceZones()[zoneI].flipMap()[index];
                }

                meshMod.modifyFace
                (
                    mesh.faces()[facei],  // modified face
                    facei,                      // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    true,                           // face flip
                    newPatchID,                         // patch for face
                    zoneI,                          // zone for face
                    zoneFlip                        // face flip in zone
                );
                nUpdated++;
            }
        }
    }

    if (returnReduce(nUpdated, sumOp<label>()) > 0)
    {
        // Change the mesh. No inflation.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }
    }

    return;
}


void Foam::autoExtrude::addZones
(
    const dictionary& edict,
    labelList& cellZoneSource,
    List<labelPair>& faceZoneSource
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    if (edict.found("cellZone"))
    {
        word cellZoneName = edict.lookup("cellZone");
        cellZoneMesh& cellZones = mesh.cellZones();
        label zoneI = cellZones.findZoneID(cellZoneName);
        if (zoneI == -1)
        {
            zoneI = cellZones.size();
            cellZones.setSize(zoneI+1);
            cellZones.set
            (
                zoneI,
                new cellZone
                (
                    cellZoneName, // name
                    labelList(0),     // addressing
                    zoneI,          // index
                    cellZones       // cellZoneMesh
                )
            );
        }
        cellZoneSource = zoneI;
    }

    if (edict.found("sourceFaceZone"))
    {
        word faceZoneName = edict.lookup("sourceFaceZone");
        faceZoneMesh& faceZones = mesh.faceZones();
        label zoneI = faceZones.findZoneID(faceZoneName);
        if (zoneI == -1)
        {
            zoneI = faceZones.size();
            faceZones.setSize(zoneI+1);
            faceZones.set
            (
                zoneI,
                new faceZone
                (
                    faceZoneName, // name
                    labelList(0),     // addressing
                    boolList(0),
                    zoneI,          // index
                    faceZones       // cellZoneMesh
                )
            );
        }
        forAll(faceZoneSource, i)
        {
            faceZoneSource[i].first() = zoneI;
        }
    }

    if (edict.found("targetFaceZone"))
    {
        word faceZoneName = edict.lookup("targetFaceZone");
        faceZoneMesh& faceZones = mesh.faceZones();
        label zoneI = faceZones.findZoneID(faceZoneName);
        if (zoneI == -1)
        {
            zoneI = faceZones.size();
            faceZones.setSize(zoneI+1);
            faceZones.set
            (
                zoneI,
                new faceZone
                (
                    faceZoneName, // name
                    labelList(0),     // addressing
                    boolList(0),
                    zoneI,          // index
                    faceZones       // cellZoneMesh
                 )
             );
        }
        forAll(faceZoneSource, i)
        {
            faceZoneSource[i].second() = zoneI;
        }
    }

    return;
}


void Foam::autoExtrude::flattenPatch
(
    const indirectPrimitivePatch& pp
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const labelList& meshPoints = pp.meshPoints();

    // regionise extruded patch
    boolList blockedFaces(pp.size(), true);
    regionSplit2D extrudeRegionIDs(mesh, pp, blockedFaces);

    const pointField& sfCtrs = pp.faceCentres();
    const vectorField& sfAreas = pp.faceAreas();
    const vectorField& sfNormals = pp.faceNormals();

    const label nExtrudeRegions = extrudeRegionIDs.nRegions();
    pointField extrudeRegionCentroids(nExtrudeRegions,vector::zero);
    pointField extrudeRegionNormals(nExtrudeRegions,vector::zero);
    scalarField extrudeRegionAreas(nExtrudeRegions,scalar(0));
    labelList pointRegions(meshPoints.size(), -1);

    forAll(extrudeRegionIDs, facei)
    {
        label regioni = extrudeRegionIDs[facei];
        scalar fArea = mag(sfAreas[facei]);
        extrudeRegionCentroids[regioni] += sfCtrs[facei]*fArea;
        extrudeRegionNormals[regioni] += sfNormals[facei]*fArea;
        extrudeRegionAreas[regioni] += fArea;
        const face& lf = pp.localFaces()[facei];
        forAll(lf,fp)
        {
            label pointi = lf[fp];
            pointRegions[pointi] = max(pointRegions[pointi], regioni);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        pointRegions,
        maxEqOp<label>(),
        label(-1)               // null value
    );

    Pstream::listCombineReduce(extrudeRegionCentroids, plusOp<vector>());
    Pstream::listCombineReduce(extrudeRegionNormals, plusOp<vector>());
    Pstream::listCombineReduce(extrudeRegionAreas, plusOp<scalar>());

    List<plane> extPlanes(nExtrudeRegions);
    forAll(extrudeRegionCentroids, regioni)
    {
        extrudeRegionCentroids[regioni]
           /= extrudeRegionAreas[regioni];
        extrudeRegionNormals[regioni]
           /= extrudeRegionAreas[regioni];
        extPlanes[regioni] = plane
        (
            extrudeRegionCentroids[regioni],
            extrudeRegionNormals[regioni]
        );
    }

    vectorField dispVec(mesh.nPoints(), vector::zero);
    forAll(meshPoints, pti)
    {
        label regioni = pointRegions[pti];
        label meshpointi = meshPoints[pti];
        const point& pt = mesh.points()[meshpointi];
        dispVec[meshpointi] = extPlanes[regioni].nearestPoint(pt) - pt;
    }

    syncTools::syncPointList
    (
        mesh,
        dispVec,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    dispVec += mesh.points();
    mesh.movePoints(dispVec);
}


void Foam::autoExtrude::extrudeSelected
(
    const dictionary& extrudeDict
)
{
    Info<< nl << "Extruding layer at selected boundary faces"<<  nl <<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    //Re-balance mesh before extruding
    preBalance(extrudeDict, false);

    forAllConstIter(dictionary, extrudeDict, iter)
    {
        const word& key = iter().keyword();

        Info<<"Extruding : "<< key <<endl;

        const dictionary& edict = extrudeDict.subDict(key);

        if (!edict.found("sourcePatches"))
        {
            WarningInFunction
                <<"No sourcePatches keyword found. "
                <<"Skipping extrusion entry : " << key
                << nl << endl;
            continue;
        }

        label nLayers = edict.lookupOrDefault<label>("nLayers", 1);

        labelHashSet sourceSet = pbm.patchSet
        (
            wordReList(edict.lookup("sourcePatches")),
            false,
            true
        );

        filterExtrudePatch(sourceSet.toc());
        autoPtr<indirectPrimitivePatch> extrudePatchPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                sourceSet.toc()
            )
        );
        const indirectPrimitivePatch& extrudePatch = extrudePatchPtr();

        if (returnReduce(extrudePatch.size(), sumOp<label>()) == 0)
        {
            WarningInFunction
                <<"Extrusion contains zero sized source patch,"
                <<" skipping entry : " << key
                << nl << endl;
            continue;
        }

        Switch flatten = edict.lookupOrDefault<Switch>("flatten", false);
        if (flatten)
        {
            flattenPatch(extrudePatch);
        }

        scalar defaultExpansion =
            edict.lookupOrDefault<scalar>("expansionRatio", 1.0);
        scalarField expansionRatio(extrudePatch.nPoints(), defaultExpansion);
        pointField displacement(extrudePatch.nPoints(), vector::zero);

        extType eType = extTypeNames.read(edict.lookup("type"));

        extAction action = extAction::ADD;
        if (edict.found("action"))
        {
            action = extActionNames.read(edict.lookup("action"));
        }

        boolList blockedSourceFaces(extrudePatch.size(), true);
        regionSplit2D sourceRegionIDs(mesh, extrudePatch, blockedSourceFaces);

        //calculate first cell displacement
        switch ( eType )
        {
            case extDirection:
            {
                calculateDirectionDisplacement
                (
                    edict,
                    extrudePatch,
                    sourceRegionIDs,
                    displacement
                );
            }
            break;

            case extNormal:
            {
                calculateNormalDisplacement
                (
                    edict,
                    extrudePatch,
                    sourceRegionIDs,
                    displacement
                );
            }
            break;

            case extAverage:
            {
                calculateAverageNormalDisplacement
                (
                    edict,
                    extrudePatch,
                    sourceRegionIDs,
                    displacement
                );
            }
            break;

            case extTarget:
            {
                if
                (
                    !calculateTargetDisplacement
                    (
                        edict,
                        extrudePatch,
                        sourceSet,
                        sourceRegionIDs,
                        nLayers,
                        displacement,
                        expansionRatio
                     )
                 )
                {
                    continue;
                }
            }
            break;
        }

        label backPatchID = -1;
        if (action != extAction::ADD)
        {
            dictionary patchInfo;
            patchInfo.set("type", wallPolyPatch::typeName);
            backPatchID = meshRefiner_.addPatch
            (
                mesh,
                edict.lookupOrDefault<word>
                (
                    "exposedPatchName",
                    key
                ),
                patchInfo
             );
        }

        labelList edgePatchID;
        labelList edgeZoneID;
        boolList edgeFlip;
        labelList inflateFaceID;
        addNewProcessorPatches
        (
            extrudePatch,
            edgePatchID,
            edgeZoneID,
            edgeFlip,
            inflateFaceID
        );

        // Topo change container.
        autoPtr<polyTopoChange> meshMod
        (
            (action == extAction::NEW)
            ? new polyTopoChange(pbm.size())
            : new polyTopoChange(mesh)
        );

        // Global face indices engine
        const globalIndex globalFaces(mesh.nFaces());

        // Determine extrudePatch.edgeFaces in global numbering (so across
        // coupled patches)
        labelListList edgeGlobalFaces
        (
            addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                extrudePatch
             )
        );

        labelList exposedPatchID(0);
        if (action != extAction::ADD)
        {
            exposedPatchID.setSize
            (
                extrudePatch.size(),
                backPatchID
            );

            polyTopoChange swapMeshMod(mesh);
            forAll(extrudePatch, i)
            {
                label meshFacei = extrudePatch.addressing()[i];

                label patchi = pbm.whichPatch(meshFacei);
                label own = mesh.faceOwner()[meshFacei];
                label nei = -1;
                if (patchi == -1)
                {
                    nei = mesh.faceNeighbour()[meshFacei];
                }

                label zoneI = mesh.faceZones().whichZone(meshFacei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    label index = mesh.faceZones()[zoneI].whichFace(meshFacei);
                    zoneFlip = mesh.faceZones()[zoneI].flipMap()[index];
                }

                swapMeshMod.modifyFace
                (
                    mesh.faces()[meshFacei].reverseFace(),  // modified face
                    meshFacei,                      // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    true,                           // face flip
                    patchi,                         // patch for face
                    zoneI,                          // zone for face
                    zoneFlip                        // face flip in zone
                );
            }

            // Change the mesh. No inflation.
            autoPtr<mapPolyMesh> map = swapMeshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
            else
            {
                // Delete mesh volumes.
                mesh.clearOut();
            }
        }

        // Layers per face
        labelList nFaceLayers(extrudePatch.size(), nLayers);
        labelList nPointLayers(extrudePatch.nPoints(), nLayers);

        //Optional create zones
        labelList cellZoneSource(extrudePatch.size(), -1);
        List<labelPair> faceZoneSource(extrudePatch.size(), labelPair(-1,-1));
        addZones(edict, cellZoneSource, faceZoneSource);

        addPatchCellLayer layerExtrude(mesh, (action == extAction::ADD));

        layerExtrude.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            expansionRatio,     // expansion ratio
            extrudePatch,       // patch faces to extrude

            edgePatchID,        // if boundary edge: patch for extruded face
            edgeZoneID,         // optional zone for extruded face
            edgeFlip,
            inflateFaceID,      // mesh face that zone/patch info is from

            exposedPatchID,     // if new mesh: patches for exposed faces
            nFaceLayers,
            nPointLayers,
            displacement,
            true,
            meshMod(),
            (eType == extTarget), //use double grading
            cellZoneSource,
            faceZoneSource
        );

        List<label> newCellLevel(0);
        List<label> newPointLevel(0);
        if (action != extAction::ADD)
        {
            const labelListList& newLayerCells = layerExtrude.layerCells();
            const labelListList& newLayerPoints = layerExtrude.addedPoints();
            const labelList& copiedPts = layerExtrude.copiedPatchPoints();
            label nAddedCells = 0;
            forAll(extrudePatch, i)
            {
                nAddedCells += newLayerCells[i].size();
            }
            label nAddedPoints = 0;
            forAll(extrudePatch.meshPoints(), i)
            {
                nAddedPoints += newLayerPoints[i].size();
            }
            label nCopiedPoints = 0;
            forAll(extrudePatch.meshPoints(), i)
            {
                if (copiedPts[i] > -1)
                {
                    nCopiedPoints++;
                }
            }

            label totalPts = nAddedPoints + nCopiedPoints;
            label totalCells = nAddedCells;

            if (action != extAction::NEW)
            {
                totalPts += mesh.nPoints();
                totalCells += mesh.nCells();
            }
            newPointLevel.setSize(totalPts,0);
            newCellLevel.setSize(totalCells,0);

            const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
            const labelList& pointLevel =
                meshRefiner_.meshCutter().pointLevel();

            if (action != extAction::NEW)
            {
                forAll(cellLevel, celli)
                {
                    newCellLevel[celli] = cellLevel[celli];
                }
                forAll(pointLevel, pointi)
                {
                    newPointLevel[pointi] = pointLevel[pointi];
                }
            }

            forAll(extrudePatch, i)
            {
                label facei = extrudePatch.addressing()[i];
                label faceLevel = meshRefiner_.meshCutter().faceLevel(facei);
                if (faceLevel < 0)
                {
                    label own = mesh.faceOwner()[facei];
                    faceLevel =  cellLevel[own];
                }
                forAll(newLayerCells[i], j)
                {
                    label newcelli = newLayerCells[i][j];
                    newCellLevel[newcelli] = faceLevel;
                }
            }

            forAll(extrudePatch.meshPoints(), i)
            {
                label pointi = extrudePatch.meshPoints()[i];
                label pLevel = pointLevel[pointi];
                forAll(newLayerPoints[i], j)
                {
                    label newpointi = newLayerPoints[i][j];
                    newPointLevel[newpointi] = pLevel;
                }
                if (copiedPts[i] > -1)
                {
                    label newpointi = copiedPts[i];
                    newPointLevel[newpointi] = pLevel;
                }
            }
        }

        // Create mesh (no inflation), return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod().changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Optionally inflate mesh
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh.setInstance(meshRefiner_.timeName());

        // Update intersection info
        meshRefiner_.updateMesh(map, labelList(0));

       if (action != extAction::ADD)
       {
           meshRefiner_.meshCutter().updateLevels(newPointLevel,newCellLevel);
       }
    }

    removeUnusedPoints();
}


void Foam::autoExtrude::cleanupFaceZones()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nFaceZones  = mesh.faceZones().size();
    if (nFaceZones > 0)
    {
        polyTopoChange zoneMeshMod(mesh);
        const PtrList<surfaceZonesInfo>& surfZones =
            meshRefiner_.surfaces().surfZones();
        //Per face zone whether free-standing
        boolList freeStanding(nFaceZones, true);
        forAll(surfZones, surfI)
        {
            word faceZoneName = surfZones[surfI].faceZoneName();
            label zoneI = mesh.faceZones().findZoneID(faceZoneName);
            if (zoneI != -1 && !surfZones[surfI].freeStanding())
            {
                freeStanding[zoneI] = false;
            }
        }

        //Already decomposed so zone faces on same processor
        label nFacesReset = 0;
        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            label faceZoneI = mesh.faceZones().whichZone(facei);
            if (faceZoneI == -1)
            {
                continue;
            }

            label own = mesh.faceOwner()[facei];
            label nei = -1;
            bool switchZone = false;

            if (patchi == -1)
            {
                nei = mesh.faceNeighbour()[facei];
                label ownZoneI = mesh.cellZones().whichZone(own);
                label neiZoneI = mesh.cellZones().whichZone(nei);
                if (ownZoneI == neiZoneI && !freeStanding[faceZoneI])
                {
                    switchZone = true;
                }
            }
            else if (!patches[patchi].coupled())
            {
                switchZone = true;
            }

            if (switchZone)
            {
                zoneMeshMod.modifyFace
                (
                    mesh.faces()[facei],  // modified face
                    facei,                // label of face
                    own,                  // owner
                    nei,                  // neighbour
                    false,                // face flip
                    patchi,               // patch for face
                    -1,                   // zone for face
                    false                 // face flip in zone
                );
                nFacesReset++;
            }
        }

        if (returnReduce(nFacesReset, sumOp<label>()) != 0)
        {
            // Change the mesh. No inflation.
            autoPtr<mapPolyMesh> zoneMap =
                zoneMeshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(zoneMap);

            // Move mesh (since morphing does not do this)
            if (zoneMap().hasMotionPoints())
            {
                mesh.movePoints(zoneMap().preMotionPoints());
            }
            else
            {
                // Delete mesh volumes.
                mesh.clearOut();
            }
        }
    }

    return;
}

void Foam::autoExtrude::extrudeZeroSizedLayer
(
    const refinementParameters& refineParams
)
{
    Info<< nl << "Extruding zero-sized layer at boundary faces"<<  nl <<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    DynamicList<label> extrudedPatches(patches.size());
    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            extrudedPatches.append(patchi);
        }
    }

    indirectPrimitivePatch extrudePatch
    (
        meshRefinement::makePatch
        (
            mesh,
            extrudedPatches
        )
    );

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges
    (
        extrudePatch.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
         )
    );

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches)
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            extrudePatch
        )
    );

    labelList edgePatchID;
    labelList edgeZoneID;
    boolList edgeFlip;
    labelList inflateFaceID;
    addNewProcessorPatches
    (
        extrudePatch,
        edgePatchID,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );

    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const shellSurfaces& shells = meshRefiner_.shells();

    const labelList maxSurfLevel = surfaces.maxLevel()
        + surfaces.proxLevelIncr();
    const label maxLevel = max(max(maxSurfLevel), shells.maxIsoLevel());

    scalar extrusionLength =  0.01* edge0Len/(1<<maxLevel);

    Info<<"Extruding a distance of : "<<extrusionLength<<endl;

    boolList excludedFaces(extrudePatch.size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        extrudePatch,
        meshEdges,
        excludedFaces,
        -0.984,
        -0.984
    );
    pointField extrudePatchPointNormals =
        eClass.calculatePointNormals(excludedFaces, 0, false);

    pointField displacement(extrudePatch.nPoints());
    forAll(displacement, pointi)
    {
        const vector& patchNormal = extrudePatchPointNormals[pointi];
        displacement[pointi] = extrusionLength*patchNormal;
    }

    // Topo change container.
    autoPtr<polyTopoChange> meshMod
    (
        new polyTopoChange(mesh)
    );

    labelList exposedPatchID(0);

    label nLayers = 1;
    // Layers per face
    labelList nFaceLayers(extrudePatch.size(), nLayers);

    // Layers per point
    labelList nPointLayers(extrudePatch.nPoints(), nLayers);
    scalarField expansionRatio(extrudePatch.nPoints(), 1.0);
    addPatchCellLayer layerExtrude(mesh, true);
    layerExtrude.setRefinement
    (
        globalFaces,
        edgeGlobalFaces,

        expansionRatio,     // expansion ratio
        extrudePatch,       // patch faces to extrude

        edgePatchID,        // if boundary edge: patch for extruded face
        edgeZoneID,         // optional zone for extruded face
        edgeFlip,
        inflateFaceID,      // mesh face that zone/patch info is from

        exposedPatchID,     // if new mesh: patches for exposed faces
        nFaceLayers,
        nPointLayers,
        displacement,
        true,
        meshMod()
    );

    //Correct cell level using face level
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    List<label> newCellLevel(cellLevel.size()+extrudePatch.size());
    forAll(cellLevel, celli)
    {
        newCellLevel[celli] = cellLevel[celli];
    }
    const labelListList& newLayerCells = layerExtrude.layerCells();
    forAll(extrudePatch, i)
    {
        label facei = extrudePatch.addressing()[i];
        label faceLevel = meshRefiner_.meshCutter().faceLevel(facei);
        if (faceLevel < 0)
        {
            label own = mesh.faceOwner()[facei];
            faceLevel =  cellLevel[own];
        }
        label newCellI = newLayerCells[i][0];
        newCellLevel[newCellI] = faceLevel;
    }

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod().changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh.setInstance(meshRefiner_.timeName());

    // Update intersection info
    meshRefiner_.updateMesh(map, labelList(0));
    meshRefiner_.meshCutter().updateCellLevels(newCellLevel);

    removeUnusedPoints();

    // Now grow the extruded cells
    optimize(refineParams);

    //Optionally use weights to balance before snapping
    scalar snapWeights = refineParams.snapWeights();
    balance(snapWeights, true);

    //Handle freestanding faceZone faces
    cleanupFaceZones();

    //Reset faceZone flipMap
    meshRefiner_.resetFaceZoneFlipMap();

    return;
}


// ************************************************************************* //
