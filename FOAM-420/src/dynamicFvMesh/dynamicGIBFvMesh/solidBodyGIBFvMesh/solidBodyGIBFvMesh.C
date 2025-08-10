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
    (c) 2015-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/solidBodyGIBFvMesh/solidBodyGIBFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyGIBFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        solidBodyGIBFvMesh,
        IOobject
    );

    addToRunTimeSelectionTable
    (
        dynamicGIBFvMesh,
        solidBodyGIBFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyGIBFvMesh::initialise()
{
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                this->time().constant(),
                "triSurface",
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    initStlPointsPtr_ = new pointField(ibMeshPtr_->points());

    if (dynamicMeshCoeffs_.found("motionFunctions"))
    {
        const dictionary& motionFDict =
            dynamicMeshCoeffs_.subDict("motionFunctions");
        isSolidBody_ = motionFDict.found("solidBodyMotionFunction");

        if (isSolidBody_)
        {
            SBMFPtr_ = solidBodyMotionFunction::New(motionFDict, this->time());
        }
        else
        {
            coorFramePtr_ = coordinateFrame::lookupNew(*this, motionFDict);
            if (!coorFramePtr_->anyDynamic())
            {
                FatalErrorInFunction
                    << "coordinateFrame " << coorFramePtr_->name()
                    << " is linked to dynamic mesh but it is not dynamic."
                    << abort(FatalError);
            }
        }
    }

    if (dynamicMeshCoeffs_.found("geometryTransformationFunctions"))
    {
        const dictionary& geoTDict =
            dynamicMeshCoeffs_.subDict("geometryTransformationFunctions");

        geoTrans_.setSize(geoTDict.size());

        label fID = 0;
        forAllConstIter(dictionary, geoTDict, iter)
        {
            if (iter().isDict())
            {
                const dictionary& subDict = iter().dict();

                geoTrans_.set
                (
                    fID,
                    geometryTransformation::New(subDict)
                );
                fID++;
            }
        }
        geoTrans_.setSize(fID);
    }

    if (dynamicMeshCoeffs_.found("solverMotionFunctions"))
    {
        const dictionary& solverFDict =
            dynamicMeshCoeffs_.subDict("solverMotionFunctions");
        motionPtr_ = motionSolver::New(*this, solverFDict);
    }

    //- move in construction to store the old points
    //- otherwise the old points will be the original points of stl
    // Inconstruction it is using solidBodyGIBFvMesh::moveSurface,
    // however later it uses virtual despatch
    tmp<pointField> newSurfPointst = moveSurface(false);
    pointField& newSurfPoints = newSurfPointst.ref();

    ibMeshPtr_->movePoints(newSurfPoints);

    this->faceZones().instance() = time().timeName();
    this->cellZones().instance() = time().timeName();

    calculateHitIndex();

    Info<< dynamicMeshCoeffs_ <<endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

tmp<vectorField> Foam::solidBodyGIBFvMesh::movePolyPatch
(
    primitivePatch& pp,
    const vectorField& interSpeed
)
{
    tmp<Field<vector>> newPpPointst
    (
        new Field<vector>
        (
            pp.points()
        )
    );
    pointField& newPpPoints = newPpPointst.ref();
    forAll(newPpPoints, pI)
    {
        newPpPoints[pI] += time().deltaTValue()*interSpeed[pI];
    }
    return newPpPointst;
}


tmp<pointField> Foam::solidBodyGIBFvMesh::moveSurface(bool updateSolver)
{
    tmp<pointField> stlPointst(new pointField(ibMeshPtr_->points()));
    pointField& stlPoints = stlPointst.ref();
    const vectorField& initStlPoints = *initStlPointsPtr_;

    if (SBMFPtr_.valid() || coorFramePtr_)
    {
        septernion transform =
            coorFramePtr_
          ? coorFramePtr_->transformation()
          : SBMFPtr_().transformation();
        if (Pstream::parRun())
        {
            List<septernion> transformList(Pstream::nProcs());
            transformList[Pstream::myProcNo()] = transform;

            // Distribute transformation
            Pstream::allGatherList(transformList);

            // Do transformation based on master processor
            transform = transformList[Pstream::masterNo()];
        }

        stlPoints = transformPoints(transform, initStlPoints);
    }

    pointField oldPoints = stlPoints;

    forAll(geoTrans_, i)
    {
        stlPoints =
            geoTrans_[i].transformPoints(oldPoints)();
    }

    if (motionPtr_.valid())
    {
        if (updateSolver)
        {
            stlPoints = motionPtr_->newPoints();
        }
        else
        {
            stlPoints = motionPtr_->curPoints();
        }
    }

    return stlPointst;
}


tmp<pointField> Foam::solidBodyGIBFvMesh::recNewPointLocation
(
    const pointField& pf0,
    const labelList& addr
)
{
    tmp<pointField> pft(new pointField(pf0.size(), vector::zero));
    pointField& pf = pft.ref();

    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    const labelList& hitInd = hitIndex();

    //- Size checks for debugging
    if (addr.size() != pf0.size())
    {
        FatalErrorInFunction
            << "indexes of the points and the pointField "
            << "are not the same size"
            << abort(FatalError);
    }

    if (this->points().size() != hitInd.size())
    {
        if (hitInd.size() != 0)
        {
            FatalErrorInFunction
                << "Inconsistent sizes between hitIndices and points"
                << abort(FatalError);
        }
        else
        {
            Info<< "hitIndex is not set." << nl
                << "triSurface constructed from a faceZone." << nl
                << "Hit index is created on the fly."
                << endl;
        }
    }

    //- we want to express the point based on the coordinates of the
    //  points of the hit triangle of the stl
    //  p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a system 3x3

    forAll(addr, pI)
    {
        label addrI = addr[pI];
        label hitIndexI = hitInd[addrI];
        const triFace& triList = triSM().triSurface::operator[](hitIndexI);

        vector p0 = pf0[pI];

        vector p_1 = surP[triList[0]];
        vector p_2 = surP[triList[1]];
        vector p_3 = surP[triList[2]];

        vector w = triList.findTriangleWeights(p0, surP0);

        pf[pI] = w[0]*p_1 + w[1]*p_2 + w[2]*p_3;
    }

    return pft;
}


void Foam::solidBodyGIBFvMesh::computeOldPositionsInUnsnappedCells
(
    pointField& recAllPoints0,
    const labelList& fscp
) const
{
    pointField currentPoints(fscp.size(), vector::zero);
    forAll(fscp, pI)
    {
        label gpI = fscp[pI];
        currentPoints[pI] = this->points()[gpI];
    }

    List<pointIndexHit> nearest;
    triSM().findNearest
    (
        currentPoints,
        scalarField(currentPoints.size(), sqr(GREAT)),
        nearest
    );
    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    forAll(currentPoints, pI)
    {
        label hitIndexI = nearest[pI].index();
        const triFace& triList = triSM().triSurface::operator[](hitIndexI);
        vector disp = surP0[triList[0]]-surP[triList[0]];

        recAllPoints0[fscp[pI]] = currentPoints[pI]+disp;
    }
}


void Foam::solidBodyGIBFvMesh::calculateHitIndex()
{
    labelIOList& hitIndex = *hitIndexPtr_;
    if (this->points().size() != hitIndex.size())
    {
        if (hitIndex.size() != 0)
        {
            FatalErrorInFunction
                << "Inconsistent sizes between hitIndeces and points"
                << abort(FatalError);
        }
        else
        {
            Info<< "hitIndex is not set."
                << "triSurface constructed from a faceZone."
                << "hit Index is created on the fly."
                << endl;

            if (hitIndex.size()==0)
            {
                hitIndex.resize(this->points().size(), -1);

                List<pointIndexHit> nearest;
                triSM().findNearest
                (
                    this->points(),
                    scalarField(this->points().size(), sqr(GREAT)),
                    nearest
                );
                forAll(nearest, pI)
                {
                    const pointIndexHit& nearI = nearest[pI];
                    hitIndex[pI] = nearI.index();
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyGIBFvMesh::solidBodyGIBFvMesh
(
    const IOobject& io,
    const word& typeN
)
:
    dynamicGIBFvMesh(io, typeN),
    initStlPointsPtr_(nullptr),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                this->time().constant(),
                "triSurface",
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    initStlPointsPtr_ = new pointField(ibMeshPtr_->points());

    this->faceZones().instance() = time().timeName();
    this->cellZones().instance() = time().timeName();

    Info<< dynamicMeshCoeffs_ <<endl;
}


Foam::solidBodyGIBFvMesh::solidBodyGIBFvMesh(const IOobject& io)
:
    dynamicGIBFvMesh(io, typeName),
    initStlPointsPtr_(nullptr),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    initialise();
}


Foam::solidBodyGIBFvMesh::solidBodyGIBFvMesh
(
    const IOobject& io,
    const dictionary dict
)
:
    dynamicGIBFvMesh(io, dict),
    initStlPointsPtr_(nullptr),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyGIBFvMesh::~solidBodyGIBFvMesh()
{
    delete initStlPointsPtr_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidBodyGIBFvMesh::updateInit(const word& zoneName)
{
    tmp<pointField> newSurfPoints = moveSurface(false);

    ibMeshPtr_->movePoints(newSurfPoints());

    dynamicGIBFvMesh::updateInit(zoneName);
}


bool Foam::solidBodyGIBFvMesh::update()
{
    storeOldTimes();
    oldPoints_ = points();
    prevPoints_ = points();
    if (coorFramePtr_)
    {
        coorFramePtr_->updateState();
    }

    tmp<pointField> newSurfPoints = moveSurface(true);

    ibMeshPtr_->movePoints(newSurfPoints());

    const fvPatch& gibPatch(this->boundary()[masterId()]);
    faceList faces = preparePatch(gibPatch);

    pointField pointsF = gibPatch.patch().localPoints();

    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    tmp<vectorField> newPpPointst = recNewPointLocation
    (
        pointsF,
        gibPatch.patch().meshPoints()
    );

    pointsF = newPpPointst();

    pp.clearGeom();

    mapGIB mapCl = mapGIB(*this, pp);

    clearOutGIBData();

//----------------------------------------------------------//

    doUpdate(mapCl, true);

    return true;
}


tmp<Foam::vectorField> Foam::solidBodyGIBFvMesh::velocityCorrect
(
    const vectorField& pc
) const
{
    FatalErrorInFunction
        << "Boundary condition of velocity "
        << "should be computed based on the "
        << "old positions of the points"
        << abort(FatalError);

    return tmp<vectorField>(new vectorField(pc.size(), vector::zero));;
}



tmp<Foam::vectorField> Foam::solidBodyGIBFvMesh::oldBoundaryLocation() const
{
    const boolList& intP = interP();
    const vectorField& cP = this->points();

    tmp<vectorField> cP0t
        (
            new vectorField(cP.size(), vector::zero)
        );
    vectorField& cP0 = cP0t.ref();

    cP0 = velocityCorrect(cP)();

    boolList oldnewintPoints = boolList(points().size(), false);
    boolList popUpPoints = boolList(points().size(), false);
    boolList marknewFaces = boolList(faces().size(), false);
    boolList markoldFaces = boolList(faces().size(), false);


    const labelList& fZone0 = fl0();
    const labelList& fZone = fl();
    forAll(fZone0, fI)
    {
        markoldFaces[fZone0[fI]] = true;
    }
    forAll(fZone, fI)
    {
        marknewFaces[fZone[fI]] = true;
    }
    forAll(this->faces(), fI)
    {
        if (marknewFaces[fI] && markoldFaces[fI])
        {
            forAll(this->faces()[fI], pI)
            {
                label pointI = this->faces()[fI][pI];
                oldnewintPoints[pointI] = true;
            }
        }
    }

    const boolList& popUpC = popUpCells();

    forAll(popUpC, cI)
    {
        if (popUpC[cI])
        {
            const labelList& cP = this->cellPoints()[cI];
            forAll(cP, pI)
            {
                popUpPoints[cP[pI]] = true;
            }
        }
    }

    boolList boundaryP(this->points().size(), false);
    forAll(this->boundary(), pI)
    {
        const polyPatch& pp = this->boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if
            (
                !isA<emptyFvPatch>(this->boundary()[pI]) &&
                !isA<wedgeFvPatch>(this->boundary()[pI])
            )
            {
                if (pp.size())
                {
                    const labelList& pPoints = pp.meshPoints();
                    forAll(pPoints, pI)
                    {
                        label gpI = pPoints[pI];
                        boundaryP[gpI] = true;
                    }
                }
            }
        }
    }

    forAll(intP, pI)
    {
        if (intP[pI])
        {
            if (popUpPoints[pI])
            {
                //- previous location if point was and is at the interface
                cP0[pI] = prevPoints_[pI];
            }
            else
            {
                //- if point was not and it is at the interface
                //- calculate old location based on the interface velocity
                cP0[pI] = cP[pI] - cP0[pI]*this->time().deltaTValue();
            }
        }
        else
        {
            //- for all the non interface points
            cP0[pI] = cP[pI];
        }
    }


    return cP0t;
}

// ************************************************************************* //
