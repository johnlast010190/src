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

#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "global/unitConversion/unitConversion.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "externalDisplacementMeshMover/medialAxisMeshMover.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "regExp.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "meshRefinement/meshRefinement.H"
#include "meshTools/meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 80;

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::layerParameters::grownUpType,
        3
    >::names[] =
    {
        "grown",
        "collapsed",
        "auto"
    };
}
const Foam::NamedEnum<Foam::layerParameters::grownUpType, 3>
    Foam::layerParameters::grownUpTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Set possible patches for collapsing layers on selected points.
// Only used by dual and extrude method
void Foam::layerParameters::setCollapsePatches
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
{
    const dictionary& layersDict = dict.subDict("layers");

    bool globalCornerCollapse
    (
        dict.lookupOrDefault<bool>("globalCornerCollapse", false)
    );

    bool globalBaffleCollapse
    (
        dict.lookupOrDefault<bool>("globalBaffleCollapse", false)
    );

    cornerCollapse_.setSize
    (
        boundaryMesh.size(),
        globalCornerCollapse
    );

    baffleCollapse_.setSize
    (
        boundaryMesh.size(),
        globalBaffleCollapse
    );

    forAll(boundaryMesh, patchI)
    {
        const word& patchName = boundaryMesh[patchI].name();
        const polyPatch& pp = boundaryMesh[patchI];

        if (!pp.coupled() && layersDict.found(patchName))
        {
            const dictionary& layerDict = layersDict.subDict(patchName);
            if (layerDict.found("baffleCollapse"))
            {
                baffleCollapse_[patchI] =
                    readBool(layerDict.lookup("baffleCollapse"));
            }
            if (layerDict.found("cornerCollapse"))
            {
                cornerCollapse_[patchI] =
                    readBool(layerDict.lookup("cornerCollapse"));
            }
        }
        else if (pp.coupled())
        {
            baffleCollapse_[patchI] = false;
            cornerCollapse_[patchI] = false;
        }
    }
}


// Set layer growth method
void Foam::layerParameters::setLayerGrowthMethod
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
{
    const dictionary& layersDict = dict.subDict("layers");
    layerSpec_.setSize(boundaryMesh.size());

    singleSidedIDs_.setSize(0);

    if (dualExtrude_)
    {
        forAll(boundaryMesh, patchI)
        {
            const word& patchName = boundaryMesh[patchI].name();
            const polyPatch& pp = boundaryMesh[patchI];

            List<word> spec(3);
            spec[0] = "numLayers";
            spec[1] = "expansionRatio";
            spec[2] = "none";

            if (!pp.coupled() && layersDict.found(patchName))
            {
                const dictionary& layerDict = layersDict.subDict(patchName);
                layerDict.readIfPresent<Switch>("reSnap",reSnap_[patchI]);
                layerDict.readIfPresent<Switch>
                (
                    "relativeSizes",
                    relativeSizes_[patchI]
                );

                layerDict.readIfPresent<Tuple2<scalar,scalar>>
                (
                    "targetExpansion",
                    targetExpansions_[patchI]
                );

                if (layerDict.found("nSurfaceLayers"))
                {
                    numLayers_[patchI] =
                        readLabel(layerDict.lookup("nSurfaceLayers"));
                    word gUp;
                    bool found = layerDict.readIfPresent("grownUp", gUp);
                    if (found)
                    {
                        if (gUp == "true" || gUp == "on")
                        {
                            numLayers_[patchI] = -1;
                            grownUp_[patchI] = GROWN;
                        }
                        else if (gUp == "false" || gUp == "off")
                        {
                            grownUp_[patchI] = COLLAPSED;
                        }
                        else if (gUp == "automatic")
                        {
                            grownUp_[patchI] = AUTO;
                        }
                    }

                    if (layerDict.found("fch") || layerDict.found("rfch"))
                    {
                        if (layerDict.found("fch"))
                        {
                            spec[1] = "fch";
                            fch_[patchI] = readScalar
                            (
                                layerDict.lookup("fch")
                            );
                        }
                        else
                        {
                            spec[1] = "rfch";
                            fch_[patchI] =
                                readScalar(layerDict.lookup("rfch"));
                        }
                    }
                    else if (layerDict.found("expansionRatio"))
                    {
                        scalar eR =
                            readScalar(layerDict.lookup("expansionRatio"));
                        if (eR < SMALL)
                        {
                            IOWarningInFunction(layersDict)
                                <<" Expansion ratio for patch "
                                << patchName
                                <<" is less than or equal to zero."
                                <<" Resetting to default value "
                                << defaultExpansionRatio_
                                << endl;
                            eR = defaultExpansionRatio_;
                        }

                        expansionRatio_[patchI] = eR;
                    }
                }
                if (layerDict.found("maxLayerThickness"))
                {
                    spec[2] = "maxLayerThickness";
                    maxLayerThickness_[patchI] = readScalar
                    (
                        layerDict.lookup("maxLayerThickness")
                    );
                }
                else if (layerDict.found("finalLayerThickness"))
                {
                    spec[2] = "finalLayerThickness";
                    finalLayerThickness_[patchI] = readScalar
                    (
                        layerDict.lookup("finalLayerThickness")
                    );
                }
            }
            layerSpec_[patchI] =  spec;
        }
    }
    else
    {
        Switch globalFixedFCH = dict.lookupOrDefault<Switch>("fixedFCH", false);
        fixedFCH_.setSize(boundaryMesh.size(), globalFixedFCH);

        label nSingleSided = 0;
        forAll(boundaryMesh, patchI)
        {
            const word& patchName = boundaryMesh[patchI].name();
            const polyPatch& pp = boundaryMesh[patchI];

            List<word> spec(3);
            spec[0] = "numLayers";
            spec[1] = "expansionRatio";
            spec[2] = "finalLayerThickness";

            if (!pp.coupled() && layersDict.found(patchName))
            {
                const dictionary& layerDict = layersDict.subDict(patchName);

                if (layerDict.found("fixedFCH"))
                {
                    fixedFCH_[patchI] = readBool(layerDict.lookup("fixedFCH"));
                }

                layerDict.readIfPresent
                (
                    "minThickness",
                    minThickness_[patchI]
                );
                if (layerDict.found("minNumLayers"))
                {
                    minLayers_[patchI] =
                        readLabel(layerDict.lookup("minNumLayers"));
                }

                layerDict.readIfPresent<Switch>("reSnap",reSnap_[patchI]);

                if (layerDict.found("nSurfaceLayers"))
                {
                    numLayers_[patchI] =
                        readLabel(layerDict.lookup("nSurfaceLayers"));
                    word gUp;
                    bool found = layerDict.readIfPresent("grownUp", gUp);
                    if (found)
                    {
                        if (gUp == "true" || gUp == "on")
                        {
                            grownUp_[patchI] = GROWN;
                        }
                        else if (gUp == "false" || gUp == "off")
                        {
                            grownUp_[patchI] = COLLAPSED;
                        }
                        else if (gUp == "automatic")
                        {
                            grownUp_[patchI] = AUTO;
                        }

                        if (grownUp_[patchI] == GROWN && numLayers_[patchI] > 0)
                        {
                            IOWarningInFunction(layersDict)
                                << "Cannot grow up following patch where  "
                                << "layer growth is set:"
                                << patchName
                                << endl;
                        }
                    }

                    spec[0] = "numLayers";
                    if (layerDict.found("expansionRatio"))
                    {
                        scalar eR =
                            readScalar(layerDict.lookup("expansionRatio"));
                        if (eR < SMALL)
                        {
                            IOWarningInFunction(layersDict)
                                <<" Expansion ratio for patch "
                                << patchName
                                <<" is less than or equal to zero."
                                <<" Resetting to default value "
                                << defaultExpansionRatio_
                                << endl;
                            eR = defaultExpansionRatio_;
                        }

                        expansionRatio_[patchI] = eR;
                        spec[1] = "expansionRatio";
                        if (layerDict.found("finalLayerThickness"))
                        {
                            spec[2] = "finalLayerThickness";
                            finalLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("finalLayerThickness")
                            );
                        }
                        else if
                        (
                            layerDict.found("fch")
                            || layerDict.found("rfch")
                        )
                        {
                            // If keyword not found try and maintain fch at wall
                            if (!dict.found("truncateFromWall"))
                            {
                                truncateFromWall_ = false;
                            }
                            if (layerDict.found("fch"))
                            {
                                spec[2] = "fch";
                                fch_[patchI] = readScalar
                                (
                                    layerDict.lookup("fch")
                                );
                            }
                            else
                            {
                                spec[2] = "rfch";
                                fch_[patchI] =
                                    readScalar(layerDict.lookup("rfch"));
                            }
                        }
                        else if (layerDict.found("maxLayerThickness"))
                        {
                            spec[2] = "maxLayerThickness";
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                            );
                        }
                    }
                    else if (layerDict.found("finalLayerThickness"))
                    {
                        finalLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("finalLayerThickness"));
                        spec[2] = "finalLayerThickness";
                        if (layerDict.found("fch") || layerDict.found("rfch"))
                        {
                            // If keyword not found try and maintain fch at wall
                            if (!dict.found("truncateFromWall"))
                            {
                                truncateFromWall_ = false;
                            }
                            if (layerDict.found("fch"))
                            {
                                spec[1] = "fch";
                                fch_[patchI] = readScalar
                                (
                                    layerDict.lookup("fch")
                                );
                            }
                            else
                            {
                                spec[1] = "rfch";
                                fch_[patchI] = readScalar
                                (
                                    layerDict.lookup("rfch")
                                );
                            }
                        }
                        else if (layerDict.found("maxLayerThickness"))
                        {
                            spec[1] = "maxLayerThickness";
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                            );
                            if
                            (
                                finalLayerThickness_[patchI]
                                > maxLayerThickness_[patchI]
                             )
                            {
                                IOWarningInFunction(layersDict)
                                    << "Final layer thickness greater than "
                                    << "maximum layer thickness for patch: "
                                    << patchName
                                    << endl;
                            }
                        }
                    }
                    else if (layerDict.found("fch") || layerDict.found("rfch"))
                    {
                        // If keyword not found try and maintain fch at wall
                        if (!dict.found("truncateFromWall"))
                        {
                            truncateFromWall_ = false;
                        }

                        if (layerDict.found("fch"))
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("fch"));
                            spec[1] = "fch";
                        }
                        else
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("rfch"));
                            spec[1] = "rfch";
                        }

                        if (layerDict.found("maxLayerThickness"))
                        {
                            spec[2] = "maxLayerThickness";
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                            );

                            if
                            (
                                fch_[patchI] > maxLayerThickness_[patchI]
                            )
                            {
                                IOWarningInFunction(layersDict)
                                    << "First cell height is greater than "
                                    << "maximum layer thickness for patch: "
                                    << patchName
                                    << endl;
                            }
                        }
                    }
                    else if (layerDict.found("maxLayerThickness"))
                    {
                        maxLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("maxLayerThickness"));
                        spec[2] = "maxLayerThickness";
                    }
                }
                else if (layerDict.found("expansionRatio"))
                {
                    //Set dummy num layers (reset later)
                    numLayers_[patchI] = 1;
                    scalar eR = readScalar(layerDict.lookup("expansionRatio"));
                    if (eR < SMALL)
                    {
                        IOWarningInFunction(layersDict)
                            <<" Expansion ratio for patch "
                            << patchName
                            <<" is less than or equal to zero."
                            <<" Resetting to default value "
                            << defaultExpansionRatio_
                            << endl;
                        eR = defaultExpansionRatio_;
                    }

                    expansionRatio_[patchI] = eR;
                    spec[1] = "expansionRatio";

                    if (layerDict.found("finalLayerThickness"))
                    {
                        finalLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("finalLayerThickness"));
                        spec[2] = "finalLayerThickness";
                        if (layerDict.found("fch") || layerDict.found("rfch"))
                        {
                            // If keyword not found try and maintain fch at wall
                            if (!dict.found("truncateFromWall"))
                            {
                                truncateFromWall_ = false;
                            }
                            if (layerDict.found("fch"))
                            {
                                spec[0] = "fch";
                                fch_[patchI] = readScalar
                                (
                                    layerDict.lookup("fch")
                                );
                            }
                            else
                            {
                                spec[0] = "rfch";
                                fch_[patchI] = readScalar
                                (
                                    layerDict.lookup("rfch")
                                );
                            }
                        }
                        else if (layerDict.found("maxLayerThickness"))
                        {
                            spec[0] = "maxLayerThickness";
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                            );
                        }
                    }
                    else if (layerDict.found("fch") || layerDict.found("rfch"))
                    {
                        // If keyword not found try and maintain fch at wall
                        if (!dict.found("truncateFromWall"))
                        {
                            truncateFromWall_ = false;
                        }
                        if (layerDict.found("fch"))
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("fch"));
                            spec[2] = "fch";
                        }
                        else
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("rfch"));
                            spec[2] = "rfch";
                        }

                        if (layerDict.found("maxLayerThickness"))
                        {
                            spec[0] = "maxLayerThickness";
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                            );
                        }
                    }
                    else if (layerDict.found("maxLayerThickness"))
                    {
                        maxLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("maxLayerThickness"));
                        spec[2] = "maxLayerThickness";
                    }
                }
                else if (layerDict.found("finalLayerThickness"))
                {
                    //Set dummy num layers (reset later)
                    numLayers_[patchI] = 1;
                    finalLayerThickness_[patchI] =
                        readScalar(layerDict.lookup("finalLayerThickness"));
                    spec[2] = "finalLayerThickness";
                    if (layerDict.found("fch") || layerDict.found("rfch"))
                    {
                        // If keyword not found try and maintain fch at wall
                        if (!dict.found("truncateFromWall"))
                        {
                            truncateFromWall_ = false;
                        }

                        if (layerDict.found("fch"))
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("fch"));
                            spec[0] = "fch";
                        }
                        else
                        {
                            fch_[patchI] = readScalar(layerDict.lookup("rfch"));
                            spec[0] = "rfch";
                        }

                        if (layerDict.found("maxLayerThickness"))
                        {
                            maxLayerThickness_[patchI] = readScalar
                            (
                                layerDict.lookup("maxLayerThickness")
                           );
                            spec[1] = "maxLayerThickness";
                        }
                    }
                    else if (layerDict.found("maxLayerThickness"))
                    {
                        maxLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("maxLayerThickness"));
                        spec[0] = "maxLayerThickness";
                        if
                        (
                            finalLayerThickness_[patchI]
                            > maxLayerThickness_[patchI]
                        )
                        {
                            IOWarningInFunction(layersDict)
                                << "Final layer thickness greater than maximum "
                                << "layer thickness for patch: "
                                << patchName
                                << endl;
                        }
                    }
                }
                else if (layerDict.found("fch") || layerDict.found("rfch"))
                {
                    // If keyword not found try and maintain fch at wall
                    if (!dict.found("truncateFromWall"))
                    {
                        truncateFromWall_ = false;
                    }
                        //Set dummy num layers (reset later)
                    numLayers_[patchI] = 1;
                    if (layerDict.found("fch"))
                    {
                        fch_[patchI] = readScalar(layerDict.lookup("fch"));
                        spec[2] = "fch";
                    }
                    else
                    {
                        fch_[patchI] = readScalar(layerDict.lookup("fch"));
                        spec[2] = "fch";
                    }

                    if (layerDict.found("maxLayerThickness"))
                    {
                        maxLayerThickness_[patchI] =
                            readScalar(layerDict.lookup("maxLayerThickness"));
                        spec[0] = "maxLayerThickness";
                    }
                }
                else if (layerDict.found("maxLayerThickness"))
                {
                    //Set dummy num layers (reset later)
                    numLayers_[patchI] = 1;
                    maxLayerThickness_[patchI] =
                        readScalar(layerDict.lookup("maxLayerThickness"));
                    spec[0] = "fch";
                    // If keyword not found try and maintain fch at wall
                    if (!dict.found("truncateFromWall"))
                    {
                        truncateFromWall_ = false;
                    }
                }
                else
                {
                    //No specification for patch so setting to collapsed
                    grownUp_[patchI] = COLLAPSED;
                }

                if (numLayers_[patchI] > 0 && layerDict.found("singleSided"))
                {
                    if (layerDict.lookup("singleSided"))
                    {
                        singleSidedIDs_.append(patchI);
                        nSingleSided++;
                    }
                }
            }

            layerSpec_[patchI] =  spec;
        }
        singleSidedIDs_.setSize(nSingleSided);
    }
}

// Establish the patches to try growing layers up
void Foam::layerParameters::setGrownUpIDs
(
    const polyBoundaryMesh& boundaryMesh
)
{
    const polyMesh& mesh = boundaryMesh.mesh();
    scalar cosPlanar = 0.9848;

    label maxBoundarySize = boundaryMesh.size();
    reduce(maxBoundarySize, maxOp<label>());

    vectorField n(maxBoundarySize, vector::zero);
    labelList nFaces(maxBoundarySize, 0);

    forAll(boundaryMesh, patchI)
    {
        const polyPatch& pp = boundaryMesh[patchI];

        if (!pp.coupled())
        {
            forAll(pp, faceI)
            {
                label meshFaceI = pp.start() + faceI;
                n[patchI] += mesh.faceAreas()[meshFaceI]
                    / (mesh.magFaceAreas()[meshFaceI] + VSMALL);
            }
            nFaces[patchI] = pp.size();
        }
    }
    Pstream::listCombineReduce(n, plusOp<vector>());
    Pstream::listCombineReduce(nFaces, plusOp<label>());

    forAll(boundaryMesh, patchI)
    {
        if (nFaces[patchI] > 0)
        {
            n[patchI] /= nFaces[patchI];
            n[patchI] /= (mag(n[patchI]) + VSMALL);
        }
    }

    scalarField nNonPlanar(maxBoundarySize, 0.);
    forAll(boundaryMesh, patchi)
    {
        const polyPatch& pp = boundaryMesh[patchi];
        if (!pp.coupled())
        {
            forAll(pp, facei)
            {
                label meshFaceI = pp.start() + facei;
                vector norm = mesh.faceAreas()[meshFaceI]
                    / (mesh.magFaceAreas()[meshFaceI] + VSMALL);
                if ((norm & n[patchi]) < cosPlanar)
                {
                    nNonPlanar[patchi] += 1.;
                }
            }
        }
    }
    Pstream::listCombineReduce(nNonPlanar, plusOp<scalar>());

    boolList planarPatch(maxBoundarySize, true);
    forAll(boundaryMesh, patchi)
    {
        if (nFaces[patchi] > 0)
        {
            if (nNonPlanar[patchi]/nFaces[patchi] > 0.01)
            {
                planarPatch[patchi] = false;
            }
        }
    }
    Pstream::listCombineReduce(planarPatch, andOp<bool>());

    List<labelList> edgePatches(mesh.nEdges());
    forAll(mesh.edges(), edgei)
    {
        const labelList& edgeFaces = mesh.edgeFaces()[edgei];

        forAll(edgeFaces, eFI)
        {
            label facei = edgeFaces[eFI];

            if (!mesh.isInternalFace(facei))
            {
                label patchi = boundaryMesh.whichPatch(facei);
                if (!boundaryMesh[patchi].coupled())
                {
                    label sz = edgePatches[edgei].size();
                    edgePatches[edgei].setSize(sz+1);
                    edgePatches[edgei][sz] = patchi;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        edgePatches,
        uniqueEqOp(),
        labelList()       // initial value
    );

    grownUpIDs_.setSize(maxBoundarySize);
    label nGrownUp = 0;

    boolList grownNbr(maxBoundarySize, false);
    boolList anyGrownNbr(maxBoundarySize, false);

    forAll(boundaryMesh, patchI)
    {
        const polyPatch& pp = boundaryMesh[patchI];

        if (!pp.coupled() && numLayers_[patchI] == 0)
        {
            const labelList& meshEdges = pp.meshEdges();

            forAll(meshEdges, edgeI)
            {
                const label meshEdgeI = meshEdges[edgeI];
                if (edgePatches[meshEdgeI].size() > 1)
                {
                    forAll(edgePatches[meshEdgeI], i)
                    {
                        const label edgePatchI = edgePatches[meshEdgeI][i];

                        if (edgePatchI != patchI)
                        {
                            if (numLayers_[edgePatchI] > 0)
                            {
                                anyGrownNbr[patchI] = true;
                            }

                            if
                            (
                                (n[patchI] & n[edgePatchI]) > cosPlanar
                                && planarPatch[edgePatchI]
                                && numLayers_[edgePatchI] > 0
                             )
                            {
                                grownNbr[patchI] = false;
                            }
                            else if
                            (
                                numLayers_[edgePatchI] > 0
                             )
                            {
                                grownNbr[patchI] = true;
                                break;
                            }
                        }
                    }
                    if (grownNbr[patchI])
                    {
                        break;
                    }
                }
            }
        }
    }

    Pstream::listCombineReduce(anyGrownNbr, orOp<bool>());

    // Normal component of normals of connected faces.
    DynamicList<label> boundaryPatches(boundaryMesh.size());

    forAll(boundaryMesh, patchI)
    {
        const polyPatch& pp = boundaryMesh[patchI];
        if (!pp.coupled())
        {
            boundaryPatches.append(patchI);
        }
    }
    boundaryPatches.shrink();

    autoPtr<indirectPrimitivePatch> ppBoundaryPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            boundaryPatches
        )
    );
    const indirectPrimitivePatch& ppBoundary = ppBoundaryPtr();

    vectorField edgeNormal(mesh.nEdges(), vector(GREAT, GREAT, GREAT));

    // Precalculate meshEdge per pp edge
    labelList ppMeshEdges(ppBoundary.nEdges());

    forAll(ppMeshEdges, patchEdgeI)
    {
        const edge& e = ppBoundary.edges()[patchEdgeI];

        label v0 = ppBoundary.meshPoints()[e[0]];
        label v1 = ppBoundary.meshPoints()[e[1]];
        ppMeshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }

    const labelListList& edgeFaces = ppBoundary.edgeFaces();
    forAll(edgeFaces, edgeI)
    {
       const labelList& eFaces = ppBoundary.edgeFaces()[edgeI];

       label meshEdgeI = ppMeshEdges[edgeI];

       forAll(eFaces, i)
       {
          nomalsCombine()
          (
             edgeNormal[meshEdgeI],
             ppBoundary.faceNormals()[eFaces[i]]
          );
       }
    }

    syncTools::syncEdgeList
    (
       mesh,
       edgeNormal,
       nomalsCombine(),
       vector(GREAT, GREAT, GREAT)  // null value
    );

    scalarField nEdgeNonConvex(maxBoundarySize, 0.);
    scalarField nEdgeConvex(maxBoundarySize, 0.);

    forAll(ppBoundary.edges(), edgeI)
    {
       label meshEdgeI = ppMeshEdges[edgeI];
       if (edgePatches[meshEdgeI].size() > 1)
       {
          forAll(edgePatches[meshEdgeI], j)
          {
             label edgePatchI = edgePatches[meshEdgeI][j];
             if (numLayers_[edgePatchI] > 0)
             {
                const labelList& eFaces = ppBoundary.edgeFaces()[edgeI];
                forAll(eFaces, faceI)
                {
                   label meshFaceI = ppBoundary.addressing()[eFaces[faceI]];
                   label patchI = mesh.boundaryMesh().whichPatch(meshFaceI);
                   if (patchI != edgePatchI && numLayers_[patchI] == 0)
                   {
                      vector norm = mesh.faceAreas()[meshFaceI]
                         / (mesh.magFaceAreas()[meshFaceI] + VSMALL);

                      const point& edgeCentre =
                          mesh.edges()[meshEdgeI].centre(mesh.points());
                      const point& faceCentre =
                          mesh.faceCentres()[meshFaceI];

                      if
                      (
                          ((norm & edgeNormal[meshEdgeI]) > cosPlanar)
                          || (((faceCentre - edgeCentre)
                               & edgeNormal[meshEdgeI]) > 0.)
                       )
                      {
                         nEdgeConvex[patchI] += 1.;
                      }
                      else
                      {
                         nEdgeNonConvex[patchI] += 1.;
                      }
                   }
                }
             }
          }
       }
    }

    Pstream::listCombineReduce(nEdgeNonConvex, plusOp<scalar>());
    Pstream::listCombineReduce(nEdgeConvex, plusOp<scalar>());

    forAll(boundaryMesh, patchI)
    {
        if (numLayers_[patchI] == 0 && nEdgeConvex[patchI] > 0)
        {
            if (nEdgeConvex[patchI] > 0.05 * nEdgeNonConvex[patchI])
            {
                grownNbr[patchI] = false;
            }
        }
    }

    Pstream::listCombineReduce(grownNbr, orOp<bool>());

    while (true)
    {
        label nUpdated = 0;
        forAll(boundaryMesh, patchI)
        {
            const polyPatch& pp = boundaryMesh[patchI];
            if (!pp.coupled() && grownNbr[patchI] && planarPatch[patchI])
            {
                const labelList& meshEdges = pp.meshEdges();

                forAll(meshEdges, edgeI)
                {
                    const label meshEdgeI = meshEdges[edgeI];
                    forAll(edgePatches[meshEdgeI], i)
                    {
                        const label edgePatchI = edgePatches[meshEdgeI][i];
                        if (edgePatchI != patchI)
                        {
                            if (numLayers_[edgePatchI] == 0)
                            {
                                if
                                (
                                    (n[patchI] & n[edgePatchI]) > cosPlanar
                                && !grownNbr[edgePatchI]
                                )
                                {
                                    grownNbr[edgePatchI] = true;
                                    nUpdated++;
                                }
                            }
                        }
                    }
                }
            }
        }

        Pstream::listCombineReduce(grownNbr, orOp<bool>());

        reduce(nUpdated, sumOp<label>());
        if (nUpdated == 0)
        {
            break;
        }
    }

    forAll(boundaryMesh, patchI)
    {
        const polyPatch& pp = boundaryMesh[patchI];

        if (grownUp_[patchI] != COLLAPSED)
        {
            if
            (
                (
                    !pp.coupled()
                    &&  planarPatch[patchI]
                    &&  grownNbr[patchI]
                    &&  numLayers_[patchI] == 0
                 )
                ||  (grownUp_[patchI] == GROWN && anyGrownNbr[patchI])
             )
            {
                grownUpIDs_[nGrownUp] = patchI;
                nGrownUp++;
                Info<<"Patch "<<pp.name()
                    <<" chosen for growing up"<<endl;
            }
            else if (!pp.coupled() && numLayers_[patchI] == 0)
            {
                if (grownUp_[patchI] == GROWN && !anyGrownNbr[patchI])
                {
                    Info<<"Patch "<<pp.name()
                        <<" cannot be grown up since has no grown neighbour"
                        <<endl;
                }
                else
                {
                    Info<<"Patch "<<pp.name()
                        <<" not chosen for growing up"
                        <<" planarPatch: "<<planarPatch[patchI]
                        <<" grownNbr: "<<grownNbr[patchI]
                        <<endl;
                }
            }
        }
    }
    grownUpIDs_.setSize(nGrownUp);

    return;
}

bool Foam::layerParameters::tryBaffleCollapse() const
{
    if (dualConcaveCollapse_ >= 190 && dualConcaveCollapse_ <= 360)
    {
        forAll(baffleCollapse_, patchi)
        {
            if (baffleCollapse_[patchi])
            {
                return true;
            }
        }
        return false;
    }
    else
    {
        return false;
    }
}


bool Foam::layerParameters::tryCornerCollapse() const
{
    forAll(cornerCollapse_, patchi)
    {
        if (cornerCollapse_[patchi])
        {
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh,
    const bool dualExtrude
)
:
    dict_(dict),
    dualExtrude_(dualExtrude),
    numLayers_(boundaryMesh.size(), 0),
    targetExpansions_
    (
        boundaryMesh.size(),
        Tuple2<scalar, scalar>(-GREAT, GREAT)
    ),
    grownUp_(boundaryMesh.size(), AUTO),
    grownUpIDs_(labelList()),
    excludedRegions_(0),
    defaultExpansionRatio_
    (
        dict.lookupOrDefault<scalar>("expansionRatio", scalar(1.25))
    ),
    expansionRatio_
    (
        boundaryMesh.size(),
        defaultExpansionRatio_ < SMALL ? scalar(1.25) : defaultExpansionRatio_
    ),
    finalLayerThickness_
    (
        boundaryMesh.size(),
        dualExtrude_ ?
        GREAT :
        dict.lookupOrDefault<scalar>("finalLayerThickness", scalar(0.4))
    ),
    minThickness_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<scalar>("minThickness", scalar(0.2))
    ),
    maxLayerThickness_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<scalar>("maxLayerThickness", GREAT)
    ),
    fch_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<scalar>("fch", 0.001)
    ),
    relativeSizes_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<Switch>("relativeSizes", true)
    ),
    reSnap_(boundaryMesh.size(), true),
    minLayers_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<label>("minNumLayers", 0)
    ),
    concaveAngle_
    (
        dict.lookupOrDefault<scalar>("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_
    (
        dict.lookupOrDefault<label>("nGrow", 0)
    ),
    nSmoothSurfaceNormals_
    (
        dict.lookupOrDefault<label>("nSmoothSurfaceNormals", 6)
    ),
    nSmoothNormals_
    (
        dict.lookupOrDefault<label>("nSmoothNormals", 3)
    ),
    nSmoothThickness_
    (
        dict.lookupOrDefault<label>("nSmoothThickness", 10)
    ),
    maxFaceThicknessRatio_
    (
        dict.lookupOrDefault<scalar>("maxFaceThicknessRatio", scalar(2.0))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos
        (
            degToRad(dict.lookupOrDefault<scalar>("minMedialAxisAngle", 90.))
         )
    ),
    nBufferCellsNoExtrude_
    (
        dict.lookupOrDefault<label>("nBufferCellsNoExtrude", 0)
    ),
    additionalReporting_
    (
        dict.lookupOrDefault<Switch>("additionalReporting", false)
    ),
    meshShrinker_
    (
        dict.lookupOrDefault<word>
        (
            "meshShrinker",
            medialAxisMeshMover::typeName
        )
    ),
    nRelaxedIter_(labelMax),
    layerRecovery_
    (
        dict.lookupOrDefault<label>("layerRecovery", 1)
    ),
    growZoneLayers_
    (
        dict.lookupOrDefault<bool>("growZoneLayers", false)
    ),
    preBalance_
    (
        dict.lookupOrDefault<bool>("rebalance", false)
    ),
    noErrors_
    (
        dict.lookupOrDefault<bool>("noErrors", false)
    ),
    splitBaffleTerminations_
    (
        dict.lookupOrDefault<bool>("splitBaffleTerminations", false)
    ),
    writeVTK_
    (
        dict.lookupOrDefault<bool>("writeVTK", false)
    ),
    projectGrownUp_
    (
        dict.lookupOrDefault<scalar>("projectGrownUp", 0.)
    ),
    maxCellDistortion_
    (
        dict.lookupOrDefault<scalar>("maxCellDistortion", 50.)
    ),
    grownUpAngleTerminateCos_
    (
        Foam::cos
        (
            degToRad(dict.lookupOrDefault<scalar>("grownUpAngleTerminate", 90.))
        )
    ),
    maxProjectionDist_
    (
        dict.lookupOrDefault<scalar>("maxProjectionDistance", GREAT)
    ),
    nSmoothPerDualLayer_
    (
        dict.lookupOrDefault<label>("nSmoothPerDualLayer", 5)
    ),
    fanAngle_
    (
        dict.lookupOrDefault<scalar>("fanAngle", 0.)
    ),
    dualConcaveCollapse_
    (
        dict.lookupOrDefault<scalar>("dualConcaveCollapse", 0.)
    ),
    layerInterfaceWeights_
    (
        dict.lookupOrDefault<scalar>("dualLayerInterfaceWeights", 0.5)
    ),
    zoneLayersScaling_
    (
        dict.lookupOrDefault<scalar>("dualZoneLayersScaling", 1.)
    ),
    extrudeBlend_
    (
        dict.lookupOrDefault<Switch>("extrudeBlend", true)
    ),
    dualOrtho_
    (
        dict.lookupOrDefault<scalar>("dualMaxOrtho", 180.0)
    ),
    dualReSnapZones_
    (
        dict.lookupOrDefault<Switch>("dualReSnapZones", true)
    ),
    fixedFCH_
    (
        boundaryMesh.size(),
        dict.lookupOrDefault<Switch>("fixedFCH", false)
    ),
    globalWarpedSplit_(dict.lookupOrDefault<scalar>("globalWarpedSplit", -1.)),
    fchWarpedSplit_(dict.lookupOrDefault<scalar>("fchWarpedSplit", -1.)),
    curvatureSplit_(dict.lookupOrDefault<Switch>("curvatureSplit", true)),
    mergeSplitCells_
    (
        dict.lookupOrDefault<Switch>("mergeOuterSplitCells", true)
    ),
    fanTetSplit_(dict.lookupOrDefault<Switch>("fanTetSplit", false)),
    maxMergePreIter_(dict.lookupOrDefault<label>("maxLayerMergePreIter", 0)),
    maxMergePostIter_(dict.lookupOrDefault<label>("maxLayerMergePostIter", 0)),
    incrementLower_(dict.lookupOrDefault<Switch>("incrementLower", false)),
    fastGeomUpdate_(dict.lookupOrDefault<Switch>("fastGeomUpdate", true)),
    fixedNormalSlip_(dict.lookupOrDefault<Switch>("fixedNormalSlip", false)),
    squishTol_(dict.lookupOrDefault<scalar>("squishTol", 0.819)),
    twoStageExtrusion_(dict.lookupOrDefault<Switch>("twoStageExtrusion", true)),
    minStretch_(dict.lookupOrDefault<scalar>("minStretch", 0.869)),
    nFCHLayers_(dict.lookupOrDefault<label>("nFCHLayers", 1))
{
    //Check patches to collapse for dual/extrude method
    setCollapsePatches(dict,boundaryMesh);

    if (dict.found("excludedRegions"))
    {
        excludedRegions_ = dict.lookup("excludedRegions");
    }

    if (layerRecovery_ <= 0)
    {
        WarningInFunction
            << "Invalid layerRecovery parameter "<<layerRecovery_
            <<" is being reset to 1."
            << endl;
        layerRecovery_ = 1;
    }

    if (nGrow_ > 0)
    {
        WarningInFunction
            << "The nGrow parameter effect has changed with respect to 1.6.x."
            << endl
            << "Please set nGrow=0 for 1.6.x behaviour."
            << endl;
    }
    if (dict.found("nRelaxedIter"))
    {
        dict.lookup("nRelaxedIter") >> nRelaxedIter_;
    }

    // Set-up optional arguments
    growUpPatches_ = dict.lookupOrDefault<bool>("growUpPatches", true);

    maxLayerIter_ = dict.lookupOrDefault<label>("maxLayerIter", 100);

    growConcaveEdge_ = dict.lookupOrDefault<bool>("growConcaveEdge", true);
    growConvexEdge_ = dict.lookupOrDefault<bool>("growConvexEdge", true);

    truncateFromWall_ = dict.lookupOrDefault<bool>("truncateFromWall", true);

    terminationStrategy_ =
        dict.lookupOrDefault<label>("terminationStrategy", 0);

    if (dict.found("featureAngleTerminate"))
    {
        scalar angle = readScalar(dict.lookup("featureAngleTerminate"));
        cosAngleTermination_ =
        (
            Foam::cos(degToRad(0.5*angle))
        );
    }
    else
    {
        cosAngleTermination_ = 0.923;
    }

    setLayerGrowthMethod
    (
        dict,
        boundaryMesh
    );

    const dictionary& layersDict = dict.subDict("layers");
    // Check whether layer specification matches any patches
    const List<keyType> wildCards = layersDict.keys(true);

    forAll(wildCards, i)
    {
        regExp re(wildCards[i]);

        bool hasMatch = false;
        forAll(boundaryMesh, patchI)
        {
            if (re.match(boundaryMesh[patchI].name()))
            {
                hasMatch = true;
                break;
            }
        }
        if (!hasMatch)
        {
            WarningInFunction
                << "Wildcard layer specification for " << wildCards[i]
                << " does not match any patch." << endl
                << "Valid patches are " << boundaryMesh.names() << endl;
        }
    }

    const List<keyType> nonWildCards = layersDict.keys(false);

    forAll(nonWildCards, i)
    {
        if (boundaryMesh.findPatchID(nonWildCards[i]) == -1)
        {
            IOWarningInFunction(layersDict)
                << "Layer specification for " << nonWildCards[i]
                << " does not match any patch." << endl
                << "Valid patches are " << boundaryMesh.names() << endl;
        }
    }
}


// ************************************************************************* //
