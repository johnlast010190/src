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
    (c) 2015-2023 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "primitives/Tuple2/Tuple2.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "distributedTriSurfaceMesh/distributedTriSurfaceMesh.H"
#include "regionSplit/regionSplit.H"
#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementParameters::refinementParameters(const dictionary& dict)
:
    maxGlobalCells_(readLabel(dict.lookup("maxGlobalCells"))),
    maxLocalCells_(readLabel(dict.lookup("maxLocalCells"))),
    minRefineCells_(readLabel(dict.lookup("minRefinementCells"))),
    nBufferLayers_(readLabel(dict.lookup("nCellsBetweenLevels"))),
    planarAngle_
    (
        dict.lookupOrDefault
        (
            "planarAngle",
            readScalar(dict.lookup("resolveFeatureAngle"))
        )
    ),
    locationsOutsideMesh_
    (
        dict.lookupOrDefault
        (
            "locationsOutsideMesh",
            pointField(0)
        )

    ),
    minZoneRegionSize_
    (
        dict.lookupOrDefault<label>
        (
            "minZoneRegionSize",
            1
        )
    ),
    faceZoneControls_(dict.subOrEmptyDict("faceZoneControls")),
    useTopologicalSnapDetection_
    (
        dict.lookupOrDefault<bool>("useTopologicalSnapDetection", true)
    ),
    rezoneProblemCells_
    (
        dict.lookupOrDefault<Switch>
        (
            "rezoneProblemCells",
            false
        )
    ),
    newZoning_
    (
        dict.lookupOrDefault<Switch>
        (
            "newZoning",
            true
        )
    ),
    mergeSingleCouples_
    (
        dict.lookupOrDefault<Switch>
        (
            "mergeSingleCouples",
            true
        )
    ),
    balanceRefine_
    (
        dict.lookupOrDefault<Switch>
        (
            "balanceThenRefine",
            true
        )
    ),
    maxLoadUnbalance_(dict.lookupOrDefault<scalar>("maxLoadUnbalance",0.1)),
    maxCellUnbalance_(dict.lookupOrDefault<label>("maxCellUnbalance", -1)),
    snapWeights_(dict.lookupOrDefault<scalar>("balanceSnapWeights",0.0)),
    nProxRays_(dict.lookupOrDefault<label>("nProxRays",6)),
    splitCells_
    (
        dict.lookupOrDefault<Switch>
        (
            "splitCells",
            false
        )
     ),
    zoneCracks_(dict.lookupOrDefault<Switch>("zoneCracks",false)),
    refineBoundaryProblemCells_
    (
        dict.lookupOrDefault<Switch>("refineBoundaryProblemCells",false)
    ),
    handleSnapProblems_
    (
        dict.lookupOrDefault<Switch>("handleSnapProblems", true)
    ),
    interfaceRefine_
    (
        dict.lookupOrDefault<Switch>("interfaceRefine", true)
    ),
    removeUnsplittable_
    (
        dict.lookupOrDefault<Switch>("removeUnsplittableCells", true)
    ),
    additionalDualRefine_
    (
        dict.lookupOrDefault<Switch>("additionalDualRefine", true)
    ),
    interZoneBaffles_
    (
        dict.lookupOrDefault<Switch>("interZoneBaffles", false)
    ),
    fullLeakChecks_
    (
        dict.lookupOrDefault<Switch>("fullLeakChecks", true)
    ),
    zoneLeakChecks_
    (
        dict.lookupOrDefault<Switch>("zoneLeakChecks", false)
    ),
    curvatureSmooth_
    (
        dict.lookupOrDefault<label>("curvatureFieldSmooth", 0)
    ),
    holeSize_
    (
        dict.lookupOrDefault<label>("closeHoleSize", 2)
    ),
    maxCornerIter_
    (
        dict.lookupOrDefault<label>("maxCornerIter", 1)
    ),
    cornerRemoveBaffles_
    (
        dict.lookupOrDefault<Switch>("cornerRemoveBaffles", true)
    ),
    removeBoundaryAnyGrownUp_
    (
        dict.lookupOrDefault<Switch>("removeBoundaryAnyGrownUp", false)
    ),
    removeBoundaryFaceZones_
    (
        dict.lookupOrDefault<Switch>("removeBoundaryFaceZone", true)
    ),
    minFeatureLength_
    (
        dict.lookupOrDefault<scalar>("minFeatureLength",-1)
    ),
    addContactCells_
    (
        dict.lookupOrDefault<Switch>("addContactCells",false)
    ),
    extrudeExtraRemoval_
    (
        dict.lookupOrDefault<Switch>("extrudeExtraRemoval",false)
    ),
    nAddedLocations_(0),
    refineOutsideGapCell_
    (
        dict.lookupOrDefault<Switch>("refineOutsideGapCell",false)
    ),
    namedLocationsRezone_
    (
        dict.lookupOrDefault<Switch>("namedLocationsRezone",false)
    ),
    mergePreExtrude_
    (
        dict.lookupOrDefault<Switch>("mergePreExtrude",false)
    ),
    filterFreeStandingHoles_
    (
        dict.lookupOrDefault<Switch>("filterFreeStandingHoles",false)
    )
{
    if (dict.found("proximityAngle"))
    {
        proximityAngle_ = mag
        (
            Foam::cos(degToRad(readScalar(dict.lookup("proximityAngle"))))
        );
    }
    else
    {
        proximityAngle_ = 0.866;
    }

    point locationInMesh;

    if (dict.readIfPresent("locationInMesh", locationInMesh))
    {
        locationsInMesh_.append(locationInMesh);
        zonesInMesh_.append("none");    // special name for no cellZone

        if (dict.found("locationsInMesh"))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot specify both specify 'locationInMesh'"
                << " and 'locationsInMesh'"
                << exit(FatalIOError);
        }
    }
    else if (dict.found("locationsInMesh"))
    {
        Istream& is(dict.lookup("locationsInMesh"));
        token nextToken(is);

        is >> nextToken;
        is >> nextToken;

        if (nextToken == token::BEGIN_LIST)
        {
            List<Tuple2<point, word>> pointsToZone;
            dict.readIfPresent("locationsInMesh", pointsToZone);
            label nZones = locationsInMesh_.size();
            locationsInMesh_.setSize(nZones+pointsToZone.size());
            zonesInMesh_.setSize(locationsInMesh_.size());

            forAll(pointsToZone, i)
            {
                locationsInMesh_[nZones] = pointsToZone[i].first();
                zonesInMesh_[nZones] = pointsToZone[i].second();
                if (zonesInMesh_[nZones] == word::null)
                {
                    zonesInMesh_[nZones] = "none";
                }
                nZones++;
            }
        }
        else
        {
            locationsInMesh_ = pointField(dict.lookup("locationsInMesh"));
            zonesInMesh_.setSize(locationsInMesh_.size(), "none");
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Cannot find either 'locationInMesh' or 'locationsInMesh' entry"
            << exit(FatalIOError);
    }

    if (dict.found("zonesToTest"))
    {
        List<Tuple2<point, word>> pointsToZone;
        dict.readIfPresent("zonesToTest", pointsToZone);
        zonesToTestLocations_.setSize(pointsToZone.size());
        zonesToTestNames_.setSize(pointsToZone.size());

        forAll(pointsToZone, i)
        {
            zonesToTestLocations_[i] = pointsToZone[i].first();
            zonesToTestNames_[i] = pointsToZone[i].second();
            if (zonesToTestNames_[i] == word::null)
            {
                zonesToTestNames_[i] = "none";
            }
        }
    }
    else
    {
        zonesToTestLocations_.setSize(0);
        zonesToTestNames_.setSize(0);
    }

    scalar featAngle(readScalar(dict.lookup("resolveFeatureAngle")));

    if (featAngle < 0 || featAngle > 180)
    {
        curvature_ = -GREAT;
    }
    else
    {
        curvature_ = Foam::cos(degToRad(featAngle));
    }

    const dictionary* subDictPtr = dict.subDictPtr
    (
        "wrapper"
    );

    if (subDictPtr)
    {
        wrap_ =
            subDictPtr->lookupOrDefault<Switch>
            (
                "wrap",
                false
             );
    }
    else
    {
        wrap_ = false;
    }

    if (dict.found("vdbDecompositionMethod"))
    {
        vdbDecompMethod_ = word(dict.lookup("vdbDecompositionMethod"));
    }
    else
    {
        vdbDecompMethod_ = "hierarchical";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::refinementParameters::isFaceZoneBaffleChecks
(
    const word& fzName
) const
{
    if (faceZoneControls_.found(fzName))
    {
        const dictionary& fzDict = faceZoneControls_.subDict(fzName);
        bool baffleChecks =
            fzDict.lookupOrDefault<bool>("baffleChecks",true);

        return baffleChecks;
    }
    else
    {
        return true;
    }
}


bool Foam::refinementParameters::isFaceZoneControlsBoundary
(
    const word& fzName
) const
{
    if (faceZoneControls_.found(fzName))
    {
        const dictionary& fzDict = faceZoneControls_.subDict(fzName);
        word faceTypeName;
        surfaceZonesInfo::faceZoneType faceType = surfaceZonesInfo::INTERNAL;
        if (fzDict.readIfPresent("faceType", faceTypeName))
        {
            faceType = surfaceZonesInfo::faceZoneTypeNames[faceTypeName];
            if (faceType == surfaceZonesInfo::BOUNDARY)
            {
                return true;
            }
        }
    }
    return false;
}


Foam::dictionary Foam::refinementParameters::getZoneInfo
(
    const word& fzName,
    surfaceZonesInfo::faceZoneType& faceType
) const
{
    dictionary patchInfo;
    patchInfo.add("type", wallPolyPatch::typeName);
    faceType = surfaceZonesInfo::INTERNAL;

    if (faceZoneControls_.found(fzName))
    {
        const dictionary& fzDict = faceZoneControls_.subDict(fzName);

        if (fzDict.found("patchInfo"))
        {
            patchInfo = fzDict.subDict("patchInfo");
        }

        word faceTypeName;
        if (fzDict.readIfPresent("faceType", faceTypeName))
        {
            faceType = surfaceZonesInfo::faceZoneTypeNames[faceTypeName];
        }
    }
    return patchInfo;
}


Foam::labelList Foam::refinementParameters::addCellZonesToMesh
(
    polyMesh& mesh
) const
{
    labelList zoneIDs(zonesInMesh_.size(), -1);
    forAll(zonesInMesh_, i)
    {
        if (zonesInMesh_[i] != word::null && zonesInMesh_[i] != "none")
        {
            zoneIDs[i] = surfaceZonesInfo::addCellZone
            (
                zonesInMesh_[i],    // name
                labelList(0),       // addressing
                mesh
            );
        }
    }
    return zoneIDs;
}


Foam::labelList Foam::refinementParameters::findCells
(
    const bool checkInsideMesh,
    const polyMesh& mesh,
    const pointField& locations
)
{
    // Force calculation of tet-diag decomposition (for use in findCell)
    (void)mesh.tetBasePtIs();

    // Global calculation engine
    globalIndex globalCells(mesh.nCells());

    // Cell label per point
    labelList cellLabels(locations.size());

    forAll(locations, i)
    {
        const point& location = locations[i];

        label localCellI = mesh.findCell(location);

        label globalCellI = -1;

        if (localCellI != -1)
        {
            globalCellI = globalCells.toGlobal(localCellI);
        }

        reduce(globalCellI, maxOp<label>());

        if (checkInsideMesh && globalCellI == -1)
        {
            FatalErrorInFunction
                << "Point " << location
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }


        label procI = globalCells.whichProcID(globalCellI);
        label procCellI = globalCells.toLocal(procI, globalCellI);

        Info<< "Found point " << location << " in cell " << procCellI
            << " on processor " << procI << endl;

        if (globalCells.isLocal(globalCellI))
        {
            cellLabels[i] = localCellI;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


Foam::labelList Foam::refinementParameters::zonedLocations
(
    const wordList& zonesInMesh
)
{
    DynamicList<label> indices(zonesInMesh.size());

    forAll(zonesInMesh, i)
    {
        if
        (
            zonesInMesh[i] != word::null
         && zonesInMesh[i] != "none"
        )
        {
            indices.append(i);
        }
    }
    return labelList(indices, true);
}


Foam::labelList Foam::refinementParameters::unzonedLocations
(
    const wordList& zonesInMesh
)
{
    DynamicList<label> indices(0);

    forAll(zonesInMesh, i)
    {
        if
        (
            zonesInMesh[i] == word::null
         || zonesInMesh[i] == "none"
        )
        {
            indices.append(i);
        }
    }
    return labelList(indices, true);
}

void Foam::refinementParameters::addToLocationsInMesh
(
    const searchableSurfaces& allGeometry,
    const dictionary& geometryDict
)
{
    DynamicList<point> addedPoints(100*allGeometry.size());

    const wordList& fileNames = allGeometry.fileNames();
    forAllConstIter(dictionary, geometryDict, iter)
    {
        const word& key = iter().keyword();
        const dictionary& dict = geometryDict.subDict(key);
        if (dict.found("transforms"))
        {
            PtrList<dictionary> transforms(dict.lookup("transforms"));
            forAll(transforms, dictI)
            {
                const dictionary& transformDict = transforms[dictI];
                const word type(transformDict.lookup("type"));
                bool addLocation = transformDict.lookupOrDefault<bool>
                (
                    "addLocation",false
                );

                if (type == "recess" && addLocation)
                {
                    label index = -1;
                    forAll(fileNames, geomi)
                    {
                        if (fileNames[geomi] == key)
                        {
                            index = geomi;
                            break;
                        }
                    }

                    if (index != -1)
                    {
                        const searchableSurface& geom = allGeometry[index];
                        if
                        (
                            isA<triSurfaceMesh>(geom)
                            && !isA<distributedTriSurfaceMesh>(geom)
                        )
                        {
                           const triSurface& ts =
                              refCast<const triSurface>(geom);
                           vector v = transformDict.lookup("direction");
                           v /= mag(v);
                           scalar t = transformDict.lookupOrDefault
                           (
                               "distance",
                               scalar(0.0)
                           );
                           boolList borderEdge(ts.nEdges(), false);
                           labelList faceZone;
                           label numZones = ts.markZones
                           (
                               borderEdge,
                               faceZone
                           );

                           vectorField zoneCentres(numZones,vector::zero);
                           scalarField zonesAreas(numZones, scalar(0));
                           forAll(faceZone, i)
                           {
                               label zonei = faceZone[i];
                               vector areaNorm = ts.faceAreas()[i];
                               if ((areaNorm & v) < 0)
                               {
                                   areaNorm = -areaNorm;
                               }
                               scalar fArea = mag(areaNorm);
                               zonesAreas[zonei] += fArea;
                               zoneCentres[zonei] +=
                                  fArea*ts.faceCentres()[i];
                           }

                           forAll(zonesAreas, zonei)
                           {
                               zoneCentres[zonei] /= zonesAreas[zonei];
                               point keepPoint = zoneCentres[zonei]
                                  - v*t;
                               addedPoints.append(keepPoint);
                           }
                        }
                    }
                }
            }
        }

        label newPts = addedPoints.size();
        if (newPts > 0)
        {
            label sz = locationsInMesh_.size();
            locationsInMesh_.setSize(sz+newPts);
            zonesInMesh_.setSize(sz+newPts);
            forAll(addedPoints, addedi)
            {
                locationsInMesh_[sz] = addedPoints[addedi];
                zonesInMesh_[sz] = word("none");
                sz++;
            }
            nAddedLocations_ += newPts;
        }
    }

    return;
}


void Foam::refinementParameters::filterLocations
(
    const fvMesh& mesh,
    const hexRef8& meshCutter,
    const labelList& globalRegion1,
    const labelList& namedSurfaceIndex
) const
{
    if (nAddedLocations_ == 0)
    {
        return;
    }

    Info<<"Filtering " << locationsInMesh_.size() <<" locations."<<endl;

    boolList blockedFace(mesh.nFaces());
    forAll(blockedFace, faceI)
    {
        if
        (
            namedSurfaceIndex.size() && namedSurfaceIndex[faceI] == -1
            && globalRegion1[faceI] == -1
        )
        {
            blockedFace[faceI] = false;
        }
        else if (namedSurfaceIndex.size() == 0 && globalRegion1[faceI] == -1)
        {
            blockedFace[faceI] = false;
        }
        else
        {
            blockedFace[faceI] = true;
        }
    }

    // No need to sync since faceToZone already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh, blockedFace);
    blockedFace.clear();

    // Force calculation of face decomposition (used in findCell)
    (void)mesh.tetBasePtIs();

    labelList keepRegions(locationsInMesh_.size(), -1);
    // For all locationsInMesh find the cell
    forAll(locationsInMesh_, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh_[i];
        label keepRegionI = -1;

        label cellI = meshRefinement::findCell
        (
            insidePoint,
            mesh,
            meshCutter
        );

        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
        }
        reduce(keepRegionI, maxOp<label>());
        keepRegions[i] = keepRegionI;
    }

    boolList setRegions(cellRegion.nRegions(), false);
    DynamicList<point> newLocations(locationsInMesh_.size());
    DynamicList<word> newZones(locationsInMesh_.size());
//    label nRemoved = 0;
    label addedStart = locationsInMesh_.size() - nAddedLocations_;
    forAll(locationsInMesh_, locationi)
    {
        label regioni = keepRegions[locationi];
        if
        (
            regioni == -1
            || (locationi > addedStart && setRegions[regioni])
        )
        {
//            nRemoved++;
            continue;
        }
        setRegions[regioni] = true;
        newLocations.append(locationsInMesh_[locationi]);
        newZones.append(zonesInMesh_[locationi]);
    }

    Info<<"Keeping " << newLocations.size() <<" filtered locations."<<endl;

    locationsInMesh_ = newLocations;
    zonesInMesh_ = newZones;
    //Reset to 0 so no further filtering
    nAddedLocations_ = 0;

    return;
}


// ************************************************************************* //
