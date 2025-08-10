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

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/layerCellHandling/layerCellDetection.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshRefinement/meshRefinement.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

namespace Foam
{

defineTypeNameAndDebug(layerCellDetection, 0);

void Foam::layerCellDetection::findCellExtremesDistances()
{
    const vectorField& faceCentres = mesh_.faceCentres();
    const vectorField& cellCentres = mesh_.cellCentres();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(mesh_.faces(), fI)
    {
        if (!isABoundaryFace(fI))
        {
            const label& ownl = own[fI];
            const label& neil = nei[fI];

            scalar ownDistance = mag(faceCentres[fI]-cellCentres[ownl]);
            scalar neiDistance = mag(faceCentres[fI]-cellCentres[neil]);

            cellMaxDistance_[ownl] = max(cellMaxDistance_[ownl], ownDistance);

            if (ownDistance < cellMinDistance_[ownl])
            {
                cellMinDistance_[ownl] = ownDistance;
                layerGrowningDirection_[ownl] =
                    faceCentres[fI]-cellCentres[ownl];
            }

            cellMaxDistance_[neil] = max(cellMaxDistance_[neil], neiDistance);

            if (neiDistance < cellMinDistance_[neil])
            {
                cellMinDistance_[neil] = neiDistance;
                layerGrowningDirection_[neil] =
                    faceCentres[fI]-cellCentres[neil];
            }

            cellMinDistance_[neil] = min(cellMinDistance_[neil], neiDistance);
        }
        else
        {
            const label& cl = own[fI];
            scalar distance = mag(faceCentres[fI]-cellCentres[cl]);
            cellMaxDistance_[cl] = max(cellMaxDistance_[cl], distance);

            if (distance < cellMinDistance_[cl])
            {
                cellMinDistance_[cl] = distance;
                layerGrowningDirection_[cl] =
                    faceCentres[fI]-cellCentres[cl];
            }
        }
    }
}

void Foam::layerCellDetection::aspectRatioCalc()
{
    findCellExtremesDistances();
    handlePentahedraCells();
    forAll(mesh_.cells(), cI)
    {
        if (cellMinDistance_[cI] > VSMALL)
        {
            aspectRatio_[cI] = cellMaxDistance_[cI]/cellMinDistance_[cI];
        }
    }
    cellMaxDistance_.clear();
    cellMinDistance_.clear();
}

void Foam::layerCellDetection::findLayerCells()
{
    if (mesh_.foundObject<volScalarField>("layerStacks"))
    {
        Info<<"Reading Layer info from mesh object"<<endl;
        const volScalarField& layerCells =
            mesh_.lookupObject<volScalarField>("layerStacks");
        forAll(mesh_.cells(), cI)
        {
            if (layerCells[cI]>-1)
            {
                isALayerCell_[cI] = true;
                layerCells_[cI] = layerCells[cI];
            }
            else
            {
                aspectRatio_[cI] = 1;
            }
        }
    }
    else
    {
        volScalarField layerCells
        (
            IOobject
            (
                "layerStacks",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        if (layerCells.empty())
        {
            Info<<"Creating layer cells manually"<<endl;
            forAll(mesh_.cells(), cI)
            {
                if (aspectRatio_[cI]>aspectRatioTrigger_)
                {
                    isALayerCell_[cI] = true;
                    layerCells_[cI] = 1;
                }
            }
            excludePyramidCells();
        }
        else
        {
            Info<<"Reading Layer info from layerStacks"<<endl;
            forAll(mesh_.cells(), cI)
            {
                if (layerCells[cI]>-1)
                {
                    isALayerCell_[cI] = true;
                    layerCells_[cI] = layerCells[cI];
                }
                else
                {
                    aspectRatio_[cI] = 1;
                }
            }
        }
    }
}

void Foam::layerCellDetection::calculateNumberOfLayers()
{
    label numberOfLayerStacks = max(layerCells_.primitiveField());
    reduce(numberOfLayerStacks, maxOp<label>());
    labelList cellsFromEachStack(numberOfLayerStacks+1, 0);
    forAll(mesh_.cells(), cI)
    {
        if (layerCells_[cI] > -1)
        {
            cellsFromEachStack[layerCells_[cI]]++;
        }
    }

    reduce(cellsFromEachStack, sumOp<labelList>());

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if (pp.coupled())
        {
            numLayers_[patchI] = 0;
        }
        else
        {
            if (pp.faceCells().size())
            {
                const label& cL = pp.faceCells()[0];
                const label& stackL = layerCells_[cL];
                if (stackL!= -1)
                {
                    numLayers_[patchI] = cellsFromEachStack[stackL];
                }
                else
                {
                    numLayers_[patchI] = 0;
                }
            }
            else
            {
                numLayers_[patchI] = 0;
            }
        }
    }

    reduce(numLayers_, maxOp<labelList>());

}

void Foam::layerCellDetection::calculateGrowingDirection()
{
    approximateGrowingDirection();
    calcSmoothFaceNormals();
    findFirstRowLayerCells();

    const labelList& numLayer = numLayers_;

    autoPtr<indirectPrimitivePatch> ppPtr = getLayerPatch();

    const indirectPrimitivePatch& pp = ppPtr();
    const labelList& meshFaces = pp.addressing();

    forAll(pp, fI)
    {
        const label& mF = meshFaces[fI];
        const label& patchI = mesh_.boundaryMesh().whichPatch(mF);
        const label& cL = mesh_.faceOwner()[mF];
        if (firstRowLayerCell_[cL])
        {
            layerStackCounter_[mF] = 1;
            layerCounter_=1;
            rowOfLayerCells_[cL] = layerStackCounter_[mF];
            if (layerCounter_ < numLayer[patchI])
            {
                buildStack(cL, patchI);
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        isALayerProcessorBoundaryFace_,
        maxEqOp<bool>()
    );

    syncTools::syncFaceList
    (
        mesh_,
        layerDirectionProcessorFaces_,
        plusEqOp<vector>()
    );

    syncTools::syncFaceList
    (
        mesh_,
        numberOflayersProcFace_,
        maxEqOp<label>()
    );

    syncTools::syncFaceList
    (
        mesh_,
        layerStackCounter_,
        maxEqOp<label>()
    );

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& polyP = mesh_.boundaryMesh()[patchI];
        if (polyP.coupled())
        {
            forAll(polyP, fI)
            {
                const label& fL = polyP.start()+fI;
                if (isALayerProcessorBoundaryFace_[fL])
                {
                    layerCounter_ = layerStackCounter_[fL];
                    const label& cL = mesh_.faceOwner()[fL];
                    if (isALayerCell_[cL] && !isAnActualLayerCell_[cL])
                    {
                        layerDirection_[cL] = layerDirectionProcessorFaces_[fL];
                        layerCounter_++;
                        rowOfLayerCells_[cL] = layerCounter_;
                        isAnActualLayerCell_[cL] = true;
                        label layerNu = numberOflayersProcFace_[fL];
                        buildStackFromProc(cL, layerNu);
                    }
                }
            }
        }
    }
}

void Foam::layerCellDetection::buildStack(const label& cI, const label& pI)
{
    const labelList& numLayer = numLayers_;

    const labelList& cellCells = mesh_.cellCells()[cI];

    scalar MAX = 0;
    label maxLabel = cI;

    forAll(cellCells, cC)
    {
        const label& cL = cellCells[cC];
        if (isALayerCell_[cL])
        {
            vector cellToCellVec =
                mesh_.cellCentres()[cL] - mesh_.cellCentres()[cI];
            cellToCellVec /= mag(cellToCellVec);
            scalar innerProduct = cellToCellVec&layerDirection_[cI];
            scalar innerProduct2 = cellToCellVec&layerDirection_[cL];
            if
            (
                innerProduct>MAX
             && innerProduct>0.5*mag(layerDirection_[cI])
             && innerProduct2>0.4
            )
            {
                MAX = innerProduct;
                maxLabel = cL;
            }
        }
    }

    if (isALayerCell_[maxLabel] && maxLabel != cI)
    {
        layerDirection_[maxLabel] = layerDirection_[cI];
        isAnActualLayerCell_[maxLabel] = true;
        layerCounter_++;
        rowOfLayerCells_[maxLabel] = layerCounter_;
        if (layerCounter_ < numLayer[pI])
        {
            buildStack(maxLabel, pI);
        }
    }

    if ((maxLabel==cI) && (layerCounter_<numLayer[pI]))
    {
        MAX = 0;
        maxLabel = -1;
        forAll(mesh_.cells()[cI], fI)
        {
            const label& fL = mesh_.cells()[cI][fI];
            vector cellToFaceCtrVec =  mesh_.faceCentres()[fL] - mesh_.cellCentres()[cI];
            cellToFaceCtrVec /= mag(cellToFaceCtrVec);
            scalar innerProduct = cellToFaceCtrVec&layerDirection_[cI];
            if (innerProduct>MAX)
            {
                MAX = innerProduct;
                maxLabel = fL;
            }
        }
        if (maxLabel!= -1)
        {
            layerDirectionProcessorFaces_[maxLabel] = layerDirection_[cI];
            isALayerProcessorBoundaryFace_[maxLabel] = true;
            layerStackCounter_[maxLabel] = rowOfLayerCells_[cI];
            numberOflayersProcFace_[maxLabel] = numLayer[pI];
        }
    }
}

void Foam::layerCellDetection::buildStackFromProc(const label& cI, const label& layerNu)
{
    const labelList& cellCells = mesh_.cellCells()[cI];

    scalar MAX = 0;
    label maxLabel = cI;
    forAll(cellCells, cC)
    {
        const label& cL = cellCells[cC];
        if (isALayerCell_[cL])
        {
            vector cellToCellVec =
                mesh_.cellCentres()[cL] - mesh_.cellCentres()[cI];
            cellToCellVec /= mag(cellToCellVec);
            scalar innerProduct = cellToCellVec&layerDirection_[cI];
            scalar innerProduct2 = cellToCellVec&layerDirection_[cL];
            if
            (
                innerProduct>MAX
             && innerProduct>0.5*mag(layerDirection_[cI])
             && innerProduct2>0.4
            )
            {
                MAX = innerProduct;
                maxLabel = cL;
            }
        }
    }

    if (isALayerCell_[maxLabel] && maxLabel != cI)
    {
        layerDirection_[maxLabel] = layerDirection_[cI];
        isAnActualLayerCell_[maxLabel] = true;
        layerCounter_++;
        rowOfLayerCells_[maxLabel] = layerCounter_;
        if (layerCounter_ < layerNu)
        {
            buildStackFromProc(maxLabel, layerNu);
        }
    }
}
void Foam::layerCellDetection::approximateGrowingDirection()
{
    const labelList& numLayer = numLayers_;

    labelHashSet pIDs;
    forAll(mesh_.boundary(), pI)
    {
        if (numLayer[pI]>0)
        {
            pIDs.insert(pI);
        }
    }

    wallDist y(mesh_, pIDs);

    layerDirection_ = fv::gaussGrad<scalar>(mesh_).calcGrad(y.y(), y.name());
}

Foam::pointField Foam::layerCellDetection::calcSmoothFaceNormals()
{
    autoPtr<indirectPrimitivePatch> ppPtr = getLayerPatch();

    const indirectPrimitivePatch& pp = ppPtr();

    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh_.edges(),
            mesh_.pointEdges()
        )
    );

    vectorField faceNormals(pp.faceNormals());
    vectorField pointNormals(pp.meshPoints().size(), vector::zero);
    //Not sure about the effect of smoothing the face normals
/*
    forAll(pp.meshPoints(), pI)
    {
        const labelList& pointFaces = pp.pointFaces()[pI];
        forAll(pointFaces, fI)
        {
            const label& fL = pointFaces[fI];
            pointNormals[pI] += faceNormals[fL];
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pp.meshPoints(),
        pointNormals,
        plusEqOp<vector>(),
        vector::zero
    );

    forAll(pp.meshPoints(), pI)
    {
        pointNormals[pI] /= mag(pointNormals[pI]);
        const labelList& pointFaces = pp.pointFaces()[pI];
        forAll(pointFaces, fI)
        {
            const label& fL = pointFaces[fI];
            faceNormals[fL] += 0.2*pointNormals[pI];
        }
    }

    forAll(faceNormals, fI)
    {
        faceNormals[fI] /= mag(faceNormals[fI]);
    }
*/
    return faceNormals;
}

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerCellDetection::getLayerPatch()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& patchToNLayers = numLayers_;

    DynamicList<label> adaptPatches(patches.size());

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && patchToNLayers[patchI] != -1)
        {
            adaptPatches.append(patchI);
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        meshRefinement::makePatch
        (
            mesh_,
            adaptPatches
        )
    );
}

void Foam::layerCellDetection::findFirstRowLayerCells()
{
    autoPtr<indirectPrimitivePatch> ppPtr = getLayerPatch();

    const indirectPrimitivePatch& pp = ppPtr();
    pointField smoothNormals = calcSmoothFaceNormals();
    const labelList& meshFaces = pp.addressing();
    forAll(pp, fI)
    {
        const label& mF = meshFaces[fI];
        const label& cL = mesh_.faceOwner()[mF];
        if (isALayerCell_[cL])
        {
            firstRowLayerCell_[cL] = true;
            isAnActualLayerCell_[cL] = true;
            layerDirection_[cL] = -smoothNormals[fI];
        }
    }
}

void Foam::layerCellDetection::handlePentahedraCells()
{
    const labelListList& cellEdges = mesh_.cellEdges();
    forAll(mesh_.cells(), cellI)
    {
        label numberOfFaces = mesh_.cells()[cellI].size();
        if (numberOfFaces==5)
        {
            cellMaxDistance_[cellI] = SMALL;
            cellMinDistance_[cellI] = GREAT;
            forAll(cellEdges[cellI], edgeI)
            {
                label edgeLabel = cellEdges[cellI][edgeI];
                const edge& edg = mesh_.edges()[edgeLabel];
                label start = edg[0];
                label end = edg[1];
                scalar edgeLength =
                    mag(mesh_.points()[start]-mesh_.points()[end]);
                cellMaxDistance_[cellI] =
                    max(cellMaxDistance_[cellI], edgeLength);
                if (edgeLength < cellMinDistance_[cellI])
                {
                    cellMinDistance_[cellI] = edgeLength;
                    layerGrowningDirection_[cellI] =
                        mesh_.points()[start]-mesh_.points()[end];
                }
            }
        }
    }
}

bool Foam::layerCellDetection::isAPyramidCell(const label& cellI)
{
    const cell& cel = mesh_.cells()[cellI];
    label numberOfFaces = cel.size();
    if (numberOfFaces == 4)
    {
        return true;
    }
    else if (numberOfFaces == 5)
    {
        label numberOfPoints = mesh_.cellPoints()[cellI].size();
        if (numberOfPoints == 5)
        {
            return true;
        }
    }
    return false;
}

void Foam::layerCellDetection::excludePyramidCells()
{
    forAll(mesh_.cells(),cellI)
    {
        if (isAPyramidCell(cellI))
        {
            excludeCell(cellI);
        }
    }
}

void Foam::layerCellDetection::excludeCell(const label& cell)
{
    isALayerCell_[cell] = false;
    layerCells_[cell] = 0;
}

bool Foam::layerCellDetection::isABoundaryFace(const label& face)
{
    if (face<mesh_.nInternalFaces())
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool Foam::layerCellDetection::isAPrismCell(const label& cellLabel)
{
    const cell& cel = mesh_.cells()[cellLabel];
    label numberOfFaces = cel.size();
    label numberOfPoints = mesh_.cellPoints()[cellLabel].size();
    if (numberOfFaces == 5 && numberOfPoints == 6)
    {
        return true;
    }
    return false;
}

void Foam::layerCellDetection::clearData()
{
    firstRowLayerCell_.clear();
    isALayerCell_.clear();
    layerGrowningDirection_.clear();
    rowOfLayerCells_.clear();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::layerCellDetection::layerCellDetection
(
    const fvMesh& mesh,
    const scalar& aspectRatioTrigger
)
:
    mesh_(mesh),
    aspectRatioTrigger_(aspectRatioTrigger),
    patches_(mesh_.boundaryMesh()),
    cellMaxDistance_(mesh_.cells().size(), SMALL),
    cellMinDistance_(mesh_.cells().size(), GREAT),
    aspectRatio_
    (
        IOobject
        (
            "aspectRatio",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("AspectRatio", dimless, 0.0),
        "zeroGradient"
    ),
    layerCells_
    (
        IOobject
        (
            "layerCells",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("layerCells", dimless, -1.)
    ),
    layerDirection_(mesh_.cells().size(), vector::zero),
    firstRowLayerCell_(mesh_.cells().size(), false),
    isALayerCell_(mesh_.cells().size(), false),
    isAnActualLayerCell_(mesh_.cells().size(), false),
    layerGrowningDirection_(mesh_.cells().size(), vector::zero),
    numLayers_(mesh_.boundaryMesh().size(), 0),
    rowOfLayerCells_(mesh_.cells().size(), 0),
    layerCounter_(0),
    isALayerProcessorBoundaryCell_(mesh_.cells().size(), false),
    isALayerProcessorBoundaryFace_(mesh_.nFaces(), false),
    layerDirectionProcessorFaces_(mesh_.faces().size(), vector::zero),
    numberOflayersProcFace_(mesh_.nFaces(), -1),
    layerStackCounter_(mesh_.nFaces(), 0)
{
    aspectRatioCalc();
    findLayerCells();
    calculateNumberOfLayers();
    calculateGrowingDirection();
    if (debug)
    {
        volScalarField layerNumbers
        (
            IOobject
            (
                "cellNumber",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("cellNumber", dimless, 0.0)
        );
        volVectorField direction
        (
            IOobject
            (
                "direction",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("cellNumber", dimless, vector::zero)
        );
        direction.primitiveFieldRef() = layerDirection_;
        layerNumbers.primitiveFieldRef() = rowOfLayerCells_;
        layerNumbers.write();
        direction.write();
        aspectRatio_.write();
        layerCells_.write();
    }
    clearData();
}


Foam::layerCellDetection::~layerCellDetection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& layerCellDetection::aspectRatio() const
{
    return aspectRatio_;
}

const bool& layerCellDetection::isALayerCell(const label& cell) const
{
    return isAnActualLayerCell_[cell];
}

const vectorField& layerCellDetection::layerDirection() const
{
    return layerDirection_;
}

// ************************************************************************* //
}//end namespace Foam
