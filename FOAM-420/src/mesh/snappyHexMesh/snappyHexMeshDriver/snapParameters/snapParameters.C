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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "snappyHexMeshDriver/snapParameters/snapParameters.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::snapParameters::calculateRegionAndFeatureSnapPatches
(
    const dictionary& dict,
    const meshRefinement& meshRefiner
)
{
    const fvMesh& mesh = meshRefiner.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    featureEdgePatches_.setSize(boundaryMesh.size(), false);
    regionSnapPatches_.setSize(boundaryMesh.size(), false);

    const labelList& meshedPatches = meshRefiner.meshedPatches();

    if (globalFE())
    {
        forAll(meshedPatches, i)
        {
            label patchI = meshedPatches[i];
            featureEdgePatches_[patchI] = true;
        }

        forAll(mesh.faceZones(), zoneI)
        {
            label mpI, spI;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner.getFaceZoneInfo
            (
                mesh.faceZones()[zoneI].name(),
                mpI,
                spI,
                fzType
            );

            if (hasInfo)
            {
                if (zoneFeatureSnapping())
                {
                    if (fzType == surfaceZonesInfo::BOUNDARY)
                    {
                        featureEdgePatches_[mpI] = true;
                        featureEdgePatches_[spI] = true;
                    }
                    else
                    {
                        featureEdgePatches_[mpI] = true;
                        featureEdgePatches_[spI] = false;
                    }
                }
                else
                {
                    featureEdgePatches_[mpI] = false;
                    featureEdgePatches_[spI] = false;
                }
            }
        }
    }

    if (globalSnapRegions())
    {
        forAll(meshedPatches, i)
        {
            label patchI = meshedPatches[i];
            regionSnapPatches_[patchI] = true;
        }

        forAll(mesh.faceZones(), zoneI)
        {
            label mpI, spI;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner.getFaceZoneInfo
            (
                mesh.faceZones()[zoneI].name(),
                mpI,
                spI,
                fzType
            );

            if (hasInfo)
            {
                if (zoneFeatureSnapping())
                {
                    if (fzType == surfaceZonesInfo::BOUNDARY)
                    {
                        regionSnapPatches_[mpI] = true;
                        regionSnapPatches_[spI] = true;
                    }
                    else
                    {
                        regionSnapPatches_[mpI] = true;
                        regionSnapPatches_[spI] = false;
                    }
                }
                else
                {
                    regionSnapPatches_[mpI] = false;
                    regionSnapPatches_[spI] = false;
                }
            }
        }
    }

    if (const dictionary* featureEdgesDictPtr = dict.subDictPtr("featureEdges"))
    {
        forAll(boundaryMesh, patchI)
        {
            const word& patchName = boundaryMesh[patchI].name();
            const polyPatch& pp = boundaryMesh[patchI];

            if (!pp.coupled() && featureEdgesDictPtr->found(patchName))
            {
                const dictionary& featureEdgeDict =
                    featureEdgesDictPtr->subDict(patchName);

                if (featureEdgeDict.found("featureSnap"))
                {
                    featureEdgePatches_[patchI] =
                        readBool(featureEdgeDict.lookup("featureSnap"));
                }

                if (featureEdgeDict.found("snapRegionBoundaries"))
                {
                    regionSnapPatches_[patchI] =
                        readBool
                        (
                            featureEdgeDict.lookup("snapRegionBoundaries")
                         );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::snapParameters::snapParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    nSmoothPatch_(readLabel(dict.lookup("nSmoothPatch"))),
    nSmoothInternal_(dict.lookupOrDefault<label>("nSmoothInternal", 0)),
    snapTol_(readScalar(dict.lookup("tolerance"))),
    nSmoothDispl_(dict.lookupOrDefault<label>("nSolveIter",200)),
    nSnap_(readLabel(dict.lookup("nRelaxIter"))),
    featureEdgePatches_(boundaryMesh.size(), false),
    regionSnapPatches_(boundaryMesh.size(), false),
    directFeatureSnapping_
    (
        dict.lookupOrDefault<bool>("directFeatureSnapping",false)
    ),
    zoneFeatureSnapping_
    (
        dict.lookupOrDefault<bool>("zoneFeatureSnapping",true)
    ),
    concaveTol_(dict.lookupOrDefault<scalar>("concaveTol",-1)),
    collapseTol_(dict.lookupOrDefault<scalar>("collapseTol",-1)),
    nPreFeatureIter_(dict.lookupOrDefault<label>("nPreFeatureIter",0)),
    nFeatureIter_(dict.lookupOrDefault<label>("nFeatureIter",50)),
    nOuterIter_(dict.lookupOrDefault<label>("nOuterIter",1)),
    globalFE_(dict.lookupOrDefault<bool>("globalFeatureEdges",true)),
    globalSnapRegions_(dict.lookupOrDefault<bool>("globalRegionSnap",false)),
    splitDegenerateCells_
    (
        dict.lookupOrDefault<bool>("splitDegenerateCells",true)
    ),
    explicitFeatureSnap_
    (
        dict.lookupOrDefault<bool>("explicitFeatureSnap", true)
    ),
    implicitFeatureSnap_
    (
        dict.lookupOrDefault<bool>("implicitFeatureSnap", false)
    ),
    multiRegionFeatureSnap_
    (
        dict.lookupOrDefault<bool>("multiRegionFeatureSnap", false)
    ),
    detectNearSurfacesSnap_
    (
        dict.lookupOrDefault<bool>("detectNearSurfacesSnap", true)
    ),
    strictRegionSnap_
    (
        dict.lookupOrDefault<bool>("strictRegionSnap", false)
    ),
    detectBaffles_(dict.lookupOrDefault<bool>("detectBaffles", true)),
    baffleFeaturePoints_
    (
        dict.lookupOrDefault<bool>("baffleFeaturePoints", false)
    ),
    releasePoints_(dict.lookupOrDefault<bool>("releasePoints", false)),
    stringFeatures_(dict.lookupOrDefault<bool>("stringFeatures", true)),
    avoidDiagonal_(dict.lookupOrDefault<bool>("avoidDiagonal", false)),
    nFaceSplitInterval_
    (
        dict.lookupOrDefault<label>("nFaceSplitInterval", labelMin)
    ),
    concaveAngle_(dict.lookupOrDefault<scalar>("concaveAngle", 45)),
    minAreaRatio_(dict.lookupOrDefault<scalar>("minAreaRatio", 0.3)),
    averageSurfaceNormal_
    (
        dict.lookupOrDefault<bool>("averageSurfaceNormal", true)
    ),
    enlargeStencil_(dict.lookupOrDefault<bool>("enlargeStencil", false)),
    smoothSnappedSurface_
    (
        dict.lookupOrDefault<bool>("smoothSnappedSurface", true)
    ),
    reorderBaffles_(dict.lookupOrDefault<bool>("reorderBaffles", false)),
    nSmoothFeatureDisp_(dict.lookupOrDefault<label>("nSmoothFeatureDisp", 1)),
    nSliverSmooths_(dict.lookupOrDefault<label>("nSliverSmooths", 10)),
    addTetsToSplitMesh_(dict.lookupOrDefault<bool>("addTetsToSplitMesh", true)),
    featureSnapChecks_(dict.lookupOrDefault<bool>("featureSnapChecks", true)),
    writeSnapVTK_(dict.lookupOrDefault<bool>("writeSnapVTK", false)),
    nAdditionalDualLayers_
    (
        dict.lookupOrDefault<label>("nAdditionalDualLayers",0)
    ),
    preZoneSnap_(dict.lookupOrDefault<bool>("preZoneSnap", false)),
    mergeBoundaryFaces_(dict.lookupOrDefault<bool>("mergeBoundaryFaces", true)),
    mergeAcrossPatches_
    (
        dict.lookupOrDefault<bool>("mergeAcrossPatches", false)
    ),
    repatchOverlapping_
    (
        dict.lookupOrDefault<bool>("repatchOverlapping", false)
    ),
    preSmoothBaffles_
    (
        dict.lookupOrDefault<bool>("preSmoothAtBaffleEdges", false)
    ),
    acuteReflexSnap_(dict.lookupOrDefault<scalar>("acuteReflexSnapAngle", 25)),
    minAcuteReflexSnap_
    (
        dict.lookupOrDefault<scalar>
        (
            "minAcuteReflexSnapAngle", 8
        )
    ),
    conformityTol_(dict.lookupOrDefault<scalar>("conformityTol", 0.2)),
    nConformitySmooths_(dict.lookupOrDefault<label>("nConformitySmooths", 10)),
    checkConformity_(dict.lookupOrDefault<bool>("checkConformity", false)),
    dualLayerRemoval_(dict.lookupOrDefault<bool>("dualLayerRemoval", false)),
    weakFeatureSnap_
    (
        dict.lookupOrDefault<Tuple2<scalar, scalar>>
        (
            "weakFeatureSnap",
            Tuple2<scalar, scalar>(-1.0, -1.0)
        )
    ),
    preMergeExtrude_(dict.lookupOrDefault<bool>("preFaceMergeExtrude", false)),
    testPerformance_(dict.lookupOrDefault<bool>("testPerformance", false))
{
    Info<< nl
        << "Snap parameters" << nl
        << "---------------" << nl
        << endl;
}


// ************************************************************************* //
