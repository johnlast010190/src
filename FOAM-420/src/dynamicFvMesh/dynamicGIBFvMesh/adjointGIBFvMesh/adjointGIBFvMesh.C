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
    (c) 2018-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/adjointGIBFvMesh/adjointGIBFvMesh.H"
#include "dynamicGIBFvMesh/deformingBodyGIBFvMesh/deformingBodyGIBFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "dynamicGIBFvMesh/movingGIBTools/surfaceAreaSmoothing/surfaceAreaSmoothing.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointGIBFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        adjointGIBFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adjointGIBFvMesh::courantScaling(scalarField& pG)
{
    const labelList& patchPoints = this->boundary()[masterGIB_].
        patch().meshPoints();
    const labelListList& pEdges =  pointEdges();

    scalarField adgLp = scalarField(patchPoints.size(), 0);

    forAll(patchPoints, pI)
    {
        label pointI = patchPoints[pI];
        const labelList& pEdgesI = pEdges[pointI];
        forAll(pEdgesI, eI)
        {
            label edgeI = pEdgesI[eI];
            const edge& edgeII = edges()[edgeI];
            label point1 = edgeII.start();
            label point2 = edgeII.end();
            scalar edgL = mag(points()[point1] - points()[point2]);
            adgLp[pI] += edgL;
        }
        adgLp[pI] /= pEdgesI.size();
    }

    scalarField pCourant = scalarField(patchPoints.size(), 0);

    scalar maxCourant =
    adjProperties().subDict("topology").
        subDict("surfaceTracking").
        lookupOrDefault<scalar>("maxCourant", 0.5);

    autoPtr<Function1<scalar>> meanCourantPtr (
        Function1<scalar>::New
        (
            "meanCourant",
            adjProperties().subDict("topology").
                subDict("surfaceTracking")
        ));
    scalar meanCourant(meanCourantPtr->value(this->time().timeOutputValue()));


    for (int i = 1; i<20; i++)
    {
        pCourant = mag(pG)/adgLp;
       // scalar meanC = gSum(mag(pCourant))/gSum(adgLp);
        scalar meanC = gSum(mag(pG))/gSum(adgLp);

        if (meanC!=0)
        {
            scalar scale = meanCourant/meanC;
            pG *= scale;

            forAll(pG, pI)
            {
                scalar maCI = mag(pG[pI])/adgLp[pI];
                if (maCI > maxCourant)
                {
                    pG[pI] = maxCourant*sign(pG[pI])*adgLp[pI];
                }
            }
            pCourant = mag(pG)/adgLp;
        }
    }

    Info<< endl;
    Info<< "Courant max: "  << gMax(pCourant) <<endl;
    Info<< "Courant mean: " << gSum(mag(pG))/gSum(adgLp) <<endl;
}



scalarField Foam::adjointGIBFvMesh::calcCurvature()
{
    scalarField curv = this->boundary()[masterGIB_].patch().pointCurvature();
    scalarField curvTmp = curv;
    const labelListList& pEdges = this->boundary()[masterGIB_].patch().pointEdges();
    const edgeList& edges = this->boundary()[masterGIB_].patch().edges();
    forAll(pEdges, pI)
    {
        const labelList& pEdge = pEdges[pI];
        forAll(pEdge, eI)
        {
            label edgeI = pEdge[eI];
            const edge& edgeII = edges[edgeI];
            label point1 = edgeII.start();
            label point2 = edgeII.end();

            if (point1!=pI)
            {
                curvTmp[pI] += curv[point1];
            }
            else
            {
                curvTmp[pI] += curv[point2];
            }
        }
        curvTmp[pI] /= pEdges[pI].size();
    }

    curv = curvTmp;

    return curv;
}


void Foam::adjointGIBFvMesh::filterBoundaryPoints
(
    scalarField& sensN,
    const indirectPolyPatch& gibPolyPatch
)
{
    const pointField& basePoints = *basePoints_;
    const label& zoneId = gibFaceZone();
    const faceZone& flist = this->faceZones()[zoneId];

    forAll(flist, fI)
    {
        label flistI = flist[fI];
        if (flistI>=nInternalFaces())
        {
            forAll(this->faces()[flistI], pI)
            {
                label gpointI = this->faces()[flistI][pI];
                label pointI = gibPolyPatch.whichPoint(gpointI);
                if (pointI!=-1)
                {
                    if
                    (
                        (this->points()[gpointI] == basePoints[gpointI]) &&
                        (sensN[pointI]>0)
                    )
                    {
                        sensN[pointI] = 0;
                    }
                }
            }
        }
    }
}


boolList Foam::adjointGIBFvMesh::findConstraintFaces()
{
    const pointField& basePoints = *basePoints_;
    const volScalarField& G = this->lookupObject<volScalarField>("G");
    const polyPatch& poly = this->boundary()[masterGIB_].patch();
    boolList conFaces(poly.size(), false);

    if (isA<indirectPolyPatch>(poly))
    {
        const indirectPolyPatch& inPoly =
            refCast<const indirectPolyPatch>(poly);

        const labelList& addr = inPoly.fAddr();
        forAll(addr, fI)
        {
            label addrI = addr[fI];
            if (addrI>=nInternalFaces())
            {
                bool moved = false;
                forAll(this->faces()[addrI], pI)
                {
                    label gpointI = this->faces()[addrI][pI];
                    if
                    (
                        (this->points()[gpointI] != basePoints[gpointI]) &&
                        (!moved)
                    )
                    {
                        moved = true;
                    }
                }
                if ((!moved) && (G.boundaryField()[masterGIB_][fI]>0))
                {
                    conFaces[fI] = true;
                }
            }
        }
    }

    return conFaces;
}


void Foam::adjointGIBFvMesh::fixConstraintPatches
(
    pointField& snapP
)
{
    const pointField& basePoints = *basePoints_;
    forAll(this->boundary(), pI)
    {
        if
        (
            this->boundary()[pI].type() == "inlet" ||
            this->boundary()[pI].type() == "outlet" ||
            this->boundary()[pI].patch().physicalType() == "inlet" ||
            this->boundary()[pI].patch().physicalType() == "outlet"
        )
        {
            const labelList& pP = this->boundary()[pI].patch().meshPoints();
            forAll(pP, pI)
            {
                snapP[pP[pI]] = basePoints[pP[pI]];
            }
        }
    }
}


void Foam::adjointGIBFvMesh::syncProcBoundaryPoints
(
    pointField& newPoints,
    const pointField& newPointsCpy
)
{
    const labelList& patchPoints = this->boundary()[masterGIB_].
        patch().meshPoints();
    pointField relDis ( newPointsCpy - newPoints );
    syncTools::syncPointList
    (
        *this,
        patchPoints,
        relDis,
        plusEqOp<vector>(),
        vector::zero
    );
    newPoints += relDis;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::adjointGIBFvMesh::computeNewPoints
(
    primitivePatch& pp,
    const scalarField& interSpeed
)
{
    PrimitivePatchInterpolation
    <
        primitivePatch
    > pInterC(pp);

    tmp<Field<scalar>> pinterSpeedtmp
    (
        new Field<scalar>
        (
            pInterC.faceToPointInterpolate(interSpeed)
        )
    );

    scalarField& pinterSpeed = pinterSpeedtmp.ref();


    const fvPatch& gibPatch(this->boundary()[masterId()]);
    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    if (false)
    {
        simpleVTKWriter a1
        (
            pp.localFaces(),
            pp.localPoints()
        );

        a1.addPointData("pinterSpeed", pinterSpeed);
    }


    autoPtr<Function1<scalar>> curvWeightPtr (
        Function1<scalar>::New
        (
            "curvatureWeight",
            adjProperties().subDict("topology").
                subDict("surfaceTracking")
        ));

    scalar curvWeight(curvWeightPtr->value(this->time().timeOutputValue()));

    scalarField curv ( curvWeight * calcCurvature() );

    tmp<Field<vector>> ppnftmp
    (
        new Field<vector>
        (
            pInterC.faceToPointInterpolate(this->boundary()[masterGIB_].nf())
        )
    );
    const vectorField& np = ppnftmp();

    scalarField pG( pinterSpeed - curv );

    filterBoundaryPoints(pG, gibPolyPatch);

    courantScaling(pG);

    vectorField pnf( pG*np );

    //--- syncing the processor boundary points ---//
    pointField pf1 = pointField(this->points().size(), vector::zero);
    pointField pf2 = pointField(this->points().size(), vector::zero);
    const labelList& patchPoints = gibPolyPatch.meshPoints();
    forAll(patchPoints, pI)
    {
        label ppI = patchPoints[pI];
        pf1[ppI] = pnf[pI];
        pf2[ppI] = pnf[pI];
    }

    applyTwoDPlanesCorrection(pf1, pf2);

    syncTools::syncPointPositions
    (
        *this,
        pf1,
        maxEqOp<point>(),           // combine op
        point(-GREAT,-GREAT,-GREAT)    // null
    );
    syncTools::syncPointPositions
    (
        *this,
        pf2,
        minEqOp<point>(),           // combine op
        point(GREAT,GREAT,GREAT)    // null
    );


    forAll(patchPoints, pI)
    {
        const label& ppI = patchPoints[pI];
        pnf[pI] = (pf1[ppI]+pf2[ppI])/2;
    }
    ///-----------------------------------------//

    tmp<Field<vector>> newPpPointst
    (
        new Field<vector>
        (
            pp.points()
            //pp.localPoints()
        )
    );

    pointField& newPpPoints = newPpPointst.ref();

    forAll(newPpPoints, pI)
    {
        newPpPoints[pI] += time().deltaTValue()*pnf[pI];
    }

    return newPpPointst;
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointGIBFvMesh::adjointGIBFvMesh(const IOobject& io)
:
    deformingBodyGIBFvMesh(io, typeName),
    adjPropertiesPtr_(nullptr)
    /*
    pGOld_
    (
        IOobject
        (
            "pGOld",
            polyMesh::instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
        nullptr
    )
        */
{
//    pGOld_ = new scalarField(this->boundary()[masterGIB_].patch().meshPoints().size(), 0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adjointGIBFvMesh::~adjointGIBFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointGIBFvMesh::updateInit()
{
    findGIBPatches();
    clearOutGIBData();
    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = this->faceZones()[zoneId];
    faceZone& fZone =  const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    DynamicList<label> dfl(this->faces().size());
    forAll(this->boundary(), pI)
    {
        if (isA<wallFvPatch>(this->boundary()[pI]))
        {
            label startPI = this->boundary()[pI].start();
            forAll(this->boundary()[pI], pfI)
            {
                dfl.append(startPI+pfI);
            }
        }
    }
    dfl.shrink();

    fZone.resize(dfl.size());
    fZone.resetAddressing(dfl, boolList(dfl.size(), false));

    updateGIB();
}


bool Foam::adjointGIBFvMesh::update()
{
    storeOldTimes();
    oldPoints_ = points();
//----------------------------------------------------------//
    DBMFPtr_->update();

    const fvPatch& gibPatch(this->boundary()[masterId()]);
    faceList faces = preparePatch(gibPatch);

    pointField pointsF = gibPatch.patch().localPoints();

    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    pointField oldp = pointsF;

    surfaceAreaSmoothing saSmoother
    (
        *this,
        adjProperties().subDict("topology").subDict("surfaceTracking"),
        masterGIB_,
        findConstraintFaces()
    );

    saSmoother.update();

    tmp<vectorField> newPpPointst =
        computeNewPoints
        (
            pp,
            saSmoother.smoothSens()()
        );

    pointsF = newPpPointst();

    nearBoundaryIntersectionsChecking
    (
        pointsF
    );

    checkConcaveBoundaryPatchPoints(pointsF);

    mapGIB mapCl = mapGIB(*this, pp, false);

    deleteDemandDrivenData(ibMeshPtr_);
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                this->time().constant(),
                "triSurface",
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mapCl.triS()
        );


//    writeInterface(gibPolyPatch);

    clearOutGIBData();
//----------------------------------------------------------//
    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = this->faceZones()[zoneId];
    faceZone& fZone =  const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    tmp<pointField> snapPt = findSnappedPoints();
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList&  fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    correctConstraintPatches(snapP);

    fixConstraintPatches(snapP);
    correctBoundaryPointsOnBaseMesh(snapP);

    syncPoints(snapP);

    fvMesh::moveGIBPoints(snapP);
    updateGIB();

    faceCellsVisDebug();

    mapCl.mapBcs();

    popShrinkFields();
    correctV0();
    popGrowFields();
    resetMeshFluxes();
    correctBCs();

    return true;
}

const IOdictionary& Foam::adjointGIBFvMesh::adjProperties()
{
    if (!adjPropertiesPtr_)
    {
        adjPropertiesPtr_ =
        &(
            this->time().lookupObject<IOdictionary>("adjointProperties")
        );
    }


    return *adjPropertiesPtr_;
}

// ************************************************************************* //
