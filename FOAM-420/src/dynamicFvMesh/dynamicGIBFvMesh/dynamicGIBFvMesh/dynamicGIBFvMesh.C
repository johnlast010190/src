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

#include "meshTools/meshTools.H"
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/volFields/volFields.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "regionSplit/regionSplit.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "cfdTools/general/include/fvCFD.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "global/unitConversion/unitConversion.H"
#include "dynamicGIBFvMesh/movingGIBTools/twoDPointGIBCorrector/twoDPointGIBCorrector.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.H"
#include "meshes/polyMesh/polyPatches/derived/inlet/inletPolyPatch.H"
#include "meshes/polyMesh/polyPatches/derived/outlet/outletPolyPatch.H"
#include "MeshedSurface/MeshedSurfaces.H"
#include "meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicGIBFvMesh, 0);
    defineRunTimeSelectionTable(dynamicGIBFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicGIBFvMesh::initialize()
{
    modifyRegion0Patch();

    findGIBPatches();
    IOobject ioBasPoints
    (
        "basePoints",
        this->time().findInstance
            (
                this->meshDir(),
                "basePoints",
                IOobject::READ_IF_PRESENT
            ),
        meshSubDir,
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    if (exists(ioBasPoints.objectPath()))
    {
        basePoints_ = new pointIOField
        (
            ioBasPoints
        );
    }
    else
    {
        basePoints_ = new pointIOField
        (
            IOobject
            (
                "basePoints",
                polyMesh::instance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        pointField& basePoints = *basePoints_;
        basePoints  = this->points();
        basePoints_->write();
    }

    IOobject ioHitIndex
    (
        "hitIndex",
        this->time().findInstance
            (
                this->meshDir(),
                "hitIndex",
                IOobject::READ_IF_PRESENT
            ),
        meshSubDir,
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );
    if (exists(ioHitIndex.objectPath()))
    {
        hitIndexPtr_ = new labelIOList
        (
            ioHitIndex
        );
    }
    else
    {
        hitIndexPtr_ = new labelIOList
        (
            IOobject
            (
                "hitIndex",
                polyMesh::instance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
    }

    baseSf_ = new vectorField(this->faces().size(), vector::zero);
    baseCf_ = new vectorField(this->faces().size(), vector::zero);
    vectorField& baseSf = *baseSf_;
    vectorField& baseCf = *baseCf_;

    scalarField magfAreas(this->faces().size());

    baseCC_ = new vectorField(this->cells().size(), vector::zero);
    vectorField& baseCC = *baseCC_;

    scalarField cellVols(scalarField(baseCC.size(), 0));
    const pointField& basePoints = *basePoints_;

    makeFaceCentresAndAreas(basePoints, baseCf, baseSf, magfAreas);
    makeCellCentresAndVols (baseCf, baseSf, baseCC, cellVols);
    makeConcavePoints();
    makeNormal2DEdges();

    this->faceZones().instance() = time().timeName();
    this->cellZones().instance() = time().timeName();
    this->faceZones().writeOpt() = IOobject::AUTO_WRITE;
    this->cellZones().writeOpt() = IOobject::AUTO_WRITE;

    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = this->faceZones()[zoneId];
    fl0Ptr_ = new labelList(cfZone);
    fm0Ptr_ = new boolList(cfZone.flipMap());

    Info<< dynamicMeshCoeffs_ <<endl;
}


void Foam::dynamicGIBFvMesh::modifyRegion0Patch()
{
    if (region0Patch_.size()==0)
    {
        Info<< "Automatic search of the patch attached to fluid: " << endl;

        label pI = 0;
        do
        {
            const polyPatch& pp = this->boundary()[pI].patch();
            if
            (
                isA<inletPolyPatch>(pp) ||
                isA<outletPolyPatch>(pp) ||
                pp.physicalType() == "inlet" ||
                pp.physicalType() == "outlet" ||
                pp.physicalType() == "opening"
            )
            {
                region0Patch_.append(pp.name());
            }
            else if (isA<processorPolyPatch>(pp))
            {
                break;
            }
            pI ++;
        } while (pI < this->boundary().size());
    }
}


void Foam::dynamicGIBFvMesh::makeFl() const
{
    if (flPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    DynamicList<label> dfl(this->faces().size());

    labelList interestFaces = checkingIntersectingFaces();
    const vectorField& baseCc = *baseCC_;
    const vectorField& baseCf = *baseCf_;

    pointField onCC = pointField(interestFaces.size(), vector::zero);
    pointField nbCC = pointField(interestFaces.size(), vector::zero);

    boolList tmpMultInterFaces(this->faces().size(), false);

    pointField neiCc(this->nFaces() - this->nInternalFaces(), vector::zero);
    pointField ownCc(this->nFaces() - this->nInternalFaces(), vector::zero);

    forAll(this->boundary(), pI)
    {
        const polyPatch& pp = this->boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            const labelUList& faceCells = pp.faceCells();
            if (pp.coupled())
            {
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    neiCc[gfI-this->nInternalFaces()] = baseCc[faceCells[pfI]];
                    ownCc[gfI-this->nInternalFaces()] = baseCc[faceCells[pfI]];
                }
            }
            else
            {
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    neiCc[gfI-this->nInternalFaces()] = baseCf[gfI];
                    ownCc[gfI-this->nInternalFaces()] = baseCc[faceCells[pfI]];
                }
            }
        }
    }

    syncTools::swapBoundaryFacePositions(*this, neiCc);

    forAll(interestFaces, fI)
    {
        label interestFacesI = interestFaces[fI];
        if (interestFacesI < this->nInternalFaces())
        {
            const label& on = this->owner()[interestFacesI];
            const label& nb = this->neighbour()[interestFacesI];
            onCC[fI] = baseCc[on];
            nbCC[fI] = baseCc[nb];
        }
        else
        {
            label pI = this->boundaryMesh().whichPatch(interestFacesI);
            const polyPatch& pp = this->boundary()[pI].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                onCC[fI] = ownCc[interestFacesI-this->nInternalFaces()];
                nbCC[fI] = neiCc[interestFacesI-this->nInternalFaces()];
            }
        }
    }

    extendBoundaryBaseMesh(interestFaces, nbCC);

    List<List<pointIndexHit>> nearest;

    if (false)
    {
        simpleVTKWriter
        (
            onCC
        ).write
        (
            "onCC_"+this->time().timeName()+".vtk"
        );
        simpleVTKWriter
        (
            nbCC
        ).write
        (
            "nbCC_"+this->time().timeName()+".vtk"
        );
    }

    triSM().findLineAll
    (
        onCC,
        nbCC,
        nearest
    );

    boolList blockedFace(this->nFaces(), false);

    forAll(interestFaces, fI)
    {
        label interestFacesI = interestFaces[fI];

        const List<pointIndexHit>& lnearI = nearest[fI];
        if (lnearI.size()>0)
        {
            if (!(lnearI.size() % 2==0))
            {
                const pointIndexHit& nearI = lnearI[0];
                if (nearI.hit())
                {
                    dfl.append(interestFacesI);
                    blockedFace[interestFacesI] = true;
                }
            }
            else if (lnearI.size() == 2)
            {
                scalar d1 = mag(onCC[fI] - nbCC[fI]);
                scalar d2 = mag(lnearI[0].hitPoint()- lnearI[1].hitPoint());
                vector n1 = triSM().faceNormals()[lnearI[0].index()];
                vector n2 = triSM().faceNormals()[lnearI[1].index()];

                bool samePoint = (lnearI[0].hitPoint()==lnearI[1].hitPoint());

                if (samePoint)
                {
                    dfl.append(interestFacesI);
                    blockedFace[interestFacesI] = true;
                }
                else
                {
                    if (d2<SMALL) //- for multiple hiting
                    {
                        if ((n1&n2) > 0)
                        {
                            dfl.append(interestFacesI);
                            blockedFace[interestFacesI] = true;
                        }
                        else
                        {
                        /*
                            Pout<< "B: Intersection for face " << fI << endl;
                            Pout<< "n1: " << n1 << tab <<
                                triSM().faceNormals()[lnearI[0].index()] << endl;
                            Pout<< "n2: " << n2 << tab <<
                                triSM().faceNormals()[lnearI[1].index()] << endl;
                                */
                        }
                    }
                    else
                    {
                    /*
                        Pout<< "C: Intersection for face " << fI << endl;
                        Pout<< "d1: " << d1 << tab << onCC[fI] << tab
                             << nbCC[fI] << endl;
                        Pout<< "d2: " << d2 << tab << lnearI[0].hitPoint()
                             << tab << lnearI[1].hitPoint() << endl;
                        Pout<< "d3: " << triSM().triSurface::operator[](lnearI[0].index()) << endl;
                        Pout<< "d4: " <<
                        triSM().triSurface::operator[](lnearI[1].index()) << endl;
                        */
                        const face& face1 = triSM().triSurface::operator[](lnearI[0].index());
                        const face& face2 = triSM().triSurface::operator[](lnearI[1].index());
                        bool tolIssue = checkIfEdgeToleranceIntersecting(face1, face2);
                        if (tolIssue && d2<(1e-10*d1))
                        {
                            dfl.append(interestFacesI);
                            blockedFace[interestFacesI] = true;
                        }
                    }
                }
                tmpMultInterFaces[interestFacesI] = true;
            }
            else
            {
                if (debugMode_)
                {
                    Pout<< "Intersection checking for face: "
                         << interestFacesI << tab
                         << "intersection info: " << endl;
                    Pout<< lnearI << endl;
                }
            }
        }
    }


    dfl.shrink();
    flPtr_ = new labelList(dfl);


    if (false)
    {
        const pointField& basePoints = *basePoints_;
        simpleVTKWriter
        (
            this->faces(),
            labelList(dfl),
            basePoints
        ).write
        (
            "fl0"+this->time().timeName()+".vtk"
        );
    }

    multInterFacesPtr_ = new boolList(tmpMultInterFaces);
}


void Foam::dynamicGIBFvMesh::makeFm() const
{
    if (fmPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    calculateFlipmap();
    if (debugMode_)
    {
        regionVisDebug();
    }
}


void Foam::dynamicGIBFvMesh::makeInterPoints() const
{
    if (interPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    interPointsPtr_ = new boolList(this->points().size(), false);
    boolList& interPoints = *interPointsPtr_;

    flipCells(interPoints, true);
}


void Foam::dynamicGIBFvMesh::makeInterPoints0() const
{
    if (interPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    interPoints0Ptr_ = new boolList(this->points().size(), false);

    boolList& interPoints0 = *interPoints0Ptr_;
    const labelList& fll = fl0();
    forAll(fll, fI)
    {
        const face& faceI = this->faces()[fll[fI]];
        forAll(faceI, pfI)
        {
            const label& pI = faceI[pfI];
            interPoints0[pI] = true;
        }
    }
}


void Foam::dynamicGIBFvMesh::makePhiB() const
{
    if (phiBPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    const pointField& cP = points();

    dimensionedScalar deltaT = this->time().deltaT();

    const vectorField& cP0 = recAllPoints0();

    tmp<scalarField> tsweptVols(new scalarField(this->faces().size()));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(this->faces(), fI)
    {
        sweptVols[fI] = this->faces()[fI].sweptVol(cP0, cP);
    }
    phiBPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "mPhi",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimVolume/dimTime
    );
    surfaceScalarField& mPhi = *phiBPtr_;

    scalarField& mPhiVol = mPhi.ref();
    surfaceScalarField::Boundary& mPhiB = mPhi.boundaryFieldRef();

    const fvPatchList& patches = boundary();


    mPhiVol = scalarField::subField(sweptVols, nInternalFaces());
    mPhiVol /= deltaT.value();

    //- needs recheck
    forAll(patches, patchI)
    {
        if (!isA<indirectPolyPatch>(patches[patchI].patch()))
        {
            mPhiB[patchI] = patches[patchI].patchSlice(sweptVols);
            mPhiB[patchI] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, fI)
    {
        const label& fII = fzNew[fI];
        if (fII<this->nInternalFaces())
        {
            mPhiB[masterGIB_][fI] = mPhi[fII];
            mPhiB[slaveGIB_][fI]  = -mPhi[fII];
            if (fmNew[fI])
            {
                mPhiB[masterGIB_][fI] *= -1;
                mPhiB[slaveGIB_][fI] *= -1;
            }
            fiBm += 1;
            fiBs += 1;
        }
        else
        {
            label patchI = this->boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(this->boundaryMesh()[patchI]))
            {
                label lpfI = fII - this->boundaryMesh()[patchI].start();
                const labelList& fcs = this->boundary()[patchI].faceCells();
                const label& fc = fcs[lpfI];

                if (cRegion()[fc] == 0)
                {
                    mPhiB[masterGIB_][fiBm] = mPhiB[patchI][lpfI];
                    fiBm += 1;
                }
                else
                {
                    mPhiB[slaveGIB_][fiBs] = mPhiB[patchI][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::makePhiTr() const
{
    if (phiTrPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "phiTr",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimVolume/dimTime
    );
    surfaceScalarField& phiTr = *phiTrPtr_;
    phiTr = phi() - phiB();
}


void Foam::dynamicGIBFvMesh::makePhiTrS() const
{
    if (phiTrSPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrSPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "phiTrS",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        *this,
        dimVolume/dimTime
    );
    surfaceScalarField& phiTrS = *phiTrSPtr_;
    scalarField& phiTrSInt = phiTrS.ref();
    surfaceScalarField::Boundary& phiTrSB = phiTrS.boundaryFieldRef();


    const boolList faceIndi0 = faceIndicator0();
    const boolList& popUpC = popUpCells();

    //- mark faces of the pop cells that they were not an interface
    //- previously. These faces have to move and get shrinked volume
    boolList ppUp = boolList(this->points().size(), false);
    forAll(this->cells(), cI)
    {
        if (popUpC[cI]==1)
        {
            const labelList& cP(this->cellPoints()[cI]);
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                ppUp[cPI] = true;
            }
            const labelList& cellFaces(this->cells()[cI]);
            forAll(cellFaces, cfI)
            {
                const label& fI = cellFaces[cfI];
                {
                    if (faceIndi0[fI]==true)
                    {
                        const labelList& faceLI = this->faces()[fI];
                        forAll(faceLI, pI)
                        {
                            ppUp[faceLI[pI]] = false;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        *this,
        ppUp,
        plusEqOp<bool>(),
        true
    );

    //- Construction of the fluxes: based on the currect interface points and
    //the old interface points we make the pointfield with all the points at the
    //interface. By doint that we have shrinked volume
    const boolList& interfacePoints = interP();
    const boolList& interfacePoints0 = interP0();

    vectorField oldPs = prevPoints_;
    vectorField newPs = prevPoints_;

    const vectorField& oldBl = recAllPoints0();


    boolList fPopIndi  = boolList(this->faces().size(), false);
    forAll(this->points(), pI)
    {
        const labelList& pFacesI = this->pointFaces()[pI];
        if (ppUp[pI])
        {
            forAll(pFacesI, fI)
            {
                const label& pFI = pFacesI[fI];
                if (pFI<this->nInternalFaces())
                {
                    fPopIndi[pFI] = true;
                }
                else
                {
                    label patchI = this->boundaryMesh().whichPatch(pFI);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                        !isA<emptyFvPatch>(this->boundary()[patchI])
                    )
                    {
                        fPopIndi[pFI] = true;
                    }
                }
            }
        }
    }

    syncTools::syncFaceList(*this, fPopIndi, orEqOp<bool>());

    forAll(fPopIndi, fI)
    {
        if (fPopIndi[fI] == true)
        {
            const labelList& fp = this->faces()[fI];
            forAll(fp, pI)
            {
                label fpI = fp[pI];
                if
                (
                    (interfacePoints0[fpI] == false)
                    &&
                    (interfacePoints[fpI] == true)
                )
                {
                    newPs[fpI] = oldBl[fpI];
                }
            }
        }
    }

    tmp<scalarField> tsweptVols(new scalarField(this->faces().size(), 0));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(fPopIndi, fI)
    {
        if (fPopIndi[fI] == true)
        {
            sweptVols[fI] = this->faces()[fI].sweptVol(oldPs, newPs);
        }
    }

    dimensionedScalar deltaT = this->time().deltaT();
    phiTrSInt = scalarField::subField(sweptVols, nInternalFaces());
    phiTrSInt /= deltaT.value();

    const fvPatchList& patches = boundary();

    //- needs recheck
    forAll(patches, patchI)
    {
        if (!isA<indirectPolyPatch>(patches[patchI].patch()))
        {
            phiTrSB[patchI] = patches[patchI].patchSlice(sweptVols);
            phiTrSB[patchI] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, fI)
    {
        const label& fII = fzNew[fI];
        if (fII<this->nInternalFaces())
        {
            phiTrSB[masterGIB_][fI] =  phiTrSInt[fII];
            phiTrSB[slaveGIB_][fI]  = -phiTrSInt[fII];
            if (fmNew[fI])
            {
                phiTrSB[masterGIB_][fI] *= -1;
                phiTrSB[slaveGIB_][fI] *= -1;
                fiBm += 1;
                fiBs += 1;
            }
        }
        else
        {
            label patchI = this->boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(this->boundaryMesh()[patchI]))
            {
                label lpfI = fII - this->boundaryMesh()[patchI].start();
                const labelList& fcs = this->boundary()[patchI].faceCells();
                const label& fc = fcs[lpfI];

                if (cRegion()[fc] == 0)
                {
                    phiTrSB[masterGIB_][fiBm] = phiTrSB[patchI][lpfI];
                    fiBm += 1;
                }
                else
                {
                    phiTrSB[slaveGIB_][fiBs] = phiTrSB[patchI][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::makePhiTrG() const
{
    if (phiTrGPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrGPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "phiTrG",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        *this,
        dimVolume/dimTime
    );
    surfaceScalarField& phiTrG = *phiTrGPtr_;
    scalarField& phiTrGInt = phiTrG.ref();
    surfaceScalarField::Boundary& phiTrGB = phiTrG.boundaryFieldRef();

    const boolList faceIndi = faceIndicator();

    const boolList& popUpC = popUpCells();

    //- mark faces of the pop cells that they were not an interface
    //- previously. These faces have to move and get shrinked volume

    boolList ppUp = boolList(this->points().size(), false);

    forAll(this->cells(), cI)
    {
        if (popUpC[cI]==1)
        {
            const labelList& cP(this->cellPoints()[cI]);
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                ppUp[cPI] = true;
            }
            const labelList& cellFaces(this->cells()[cI]);
            forAll(cellFaces, cfI)
            {
                const label& fI = cellFaces[cfI];
                if (faceIndi[fI]==true)
                {
                    const labelList& faceLI = this->faces()[fI];
                    forAll(faceLI, pI)
                    {
                        ppUp[faceLI[pI]] = false;
                    }
                }
            }
        }
    }


    syncTools::syncPointList
    (
        *this,
        ppUp,
        plusEqOp<bool>(),
        true
    );


    boolList fPopIndi  = boolList(this->faces().size(), false);
    forAll(this->points(), pI)
    {
        const labelList& pFacesI = this->pointFaces()[pI];
        if (ppUp[pI])
        {
            forAll(pFacesI, fI)
            {
                const label& pFI = pFacesI[fI];
                if (pFI<this->nInternalFaces())
                {
                    fPopIndi[pFI] = true;
                }
                else
                {
                    label patchI = this->boundaryMesh().whichPatch(pFI);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                        !isA<emptyFvPatch>(this->boundary()[patchI])
                    )
                    {
                        fPopIndi[pFI] = true;
                    }
                }
            }
        }
    }


    syncTools::syncFaceList(*this, fPopIndi, orEqOp<bool>());

    //- Construction of the fluxes: based on the currect interface points and
    //the old interface points we make the pointfield with all the points at the
    //interface. By doint that we have shrinked volume
    //const boolList& interfacePoints = interP();
    const boolList& interfacePoints0 = interP0();

    vectorField oldPs = recAllPoints0();
    vectorField newPs = recAllPoints0();

    forAll(fPopIndi, fI)
    {
        if (fPopIndi[fI] == true)
        {
            const labelList& fp = this->faces()[fI];
            forAll(fp, pI)
            {
                label fpI = fp[pI];
                if (interfacePoints0[fpI] == true)
                {
                    newPs[fpI] = prevPoints_[fpI];
                }
            }
        }
    }

    tmp<scalarField> tsweptVols(new scalarField(this->faces().size(), 0));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(fPopIndi, fI)
    {
        if (fPopIndi[fI] == true)
        {
            sweptVols[fI] = this->faces()[fI].sweptVol(newPs, oldPs);
        }
    }

    dimensionedScalar deltaT = this->time().deltaT();
    phiTrGInt = scalarField::subField(sweptVols, nInternalFaces());
    phiTrGInt /= deltaT.value();

    const fvPatchList& patches = boundary();

    //- needs recheck
    forAll(patches, patchI)
    {
        if (!isA<indirectPolyPatch>(patches[patchI].patch()))
        {
            phiTrGB[patchI] = patches[patchI].patchSlice(sweptVols);
            phiTrGB[patchI] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, fI)
    {
        const label& fII = fzNew[fI];
        if (fII<this->nInternalFaces())
        {
            phiTrGB[masterGIB_][fI] =  phiTrGInt[fII];
            phiTrGB[slaveGIB_][fI]  = -phiTrGInt[fII];
            if (fmNew[fI])
            {
                phiTrGB[masterGIB_][fI] *= -1;
                phiTrGB[slaveGIB_][fI] *= -1;
            }
            fiBm += 1;
            fiBs += 1;
        }
        else
        {
            label patchI = this->boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(this->boundaryMesh()[patchI]))
            {
                label lpfI = fII - this->boundaryMesh()[patchI].start();
                const labelList& fcs = this->boundary()[patchI].faceCells();
                const label& fc = fcs[lpfI];

                if (cRegion()[fc] == 0)
                {
                    phiTrGB[masterGIB_][fiBm] = phiTrGB[patchI][lpfI];
                    fiBm += 1;
                }
                else
                {
                    phiTrGB[slaveGIB_][fiBs] = phiTrGB[patchI][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
//    volScalarField phiG("phiG", fvc::div(phiTrG));
//    phiG.write();
}


void Foam::dynamicGIBFvMesh::makeFaceIndicator() const
{
    if (faceIndicatorPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    faceIndicatorPtr_ = new boolList(this->faces().size(), false);
    boolList& faceIndicator = *faceIndicatorPtr_;

    const labelList& fZone = fl();

    //- treatment of the pop up/in cells
    forAll(fZone, fI)
    {
        const label& fzI = fZone[fI];
        faceIndicator[fzI] =  true;
        if (!(fzI<this->nInternalFaces()))
        {
            label patchI = this->boundaryMesh().whichPatch(fzI);
            if
            (
                isA<wedgeFvPatch>(this->boundary()[patchI]) ||
                isA<emptyFvPatch>(this->boundary()[patchI])
            )
            {
                faceIndicator[fzI] = false;
            }
        }
    }

    syncTools::syncFaceList(*this, faceIndicator, orEqOp<bool>());
}


void Foam::dynamicGIBFvMesh::makeFaceIndicator0() const
{
    if (faceIndicator0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    faceIndicator0Ptr_ = new boolList(this->faces().size(), false);
    boolList& faceIndicator0 = *faceIndicator0Ptr_;

    const labelList& fZone0 = fl0();
    forAll(fZone0, fI)
    {
        const label& fzI = fZone0[fI];
        faceIndicator0[fzI] =  true;
        if (!(fzI<this->nInternalFaces()))
        {
            label patchI = this->boundaryMesh().whichPatch(fzI);
            if
            (
                !this->boundary()[patchI].coupled()
            )
            {
                faceIndicator0[fzI] = false;
            }
        }
    }

    syncTools::syncFaceList(*this, faceIndicator0, orEqOp<bool>());
}



void Foam::dynamicGIBFvMesh::makecRegion() const
{
    if (cRegionPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList blockedFace(this->nFaces(), false);

    const labelList& fList = fl();
    forAll(fList, fI)
    {
        const label& fII = fList[fI];
        blockedFace[fII] = true;
    }

    //- parallel sync !imporant for regionSplit
    syncTools::syncFaceList(*this, blockedFace, orEqOp<bool>());

    regionSplit cellMark(*this, blockedFace);
    labelList cellIndi(cellMark);
    modifyRegionLabels(cellIndi);
    cRegionPtr_ = new labelList(cellIndi);
}

void Foam::dynamicGIBFvMesh::makecRegion0() const
{
    if (cRegion0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList blockedFace(this->nFaces(), false);

    const labelList& fList = fl0();
    forAll(fList, fI)
    {
        const label& fII = fList[fI];
        blockedFace[fII] = true;
    }

    //- parallel sync !imporant for regionSplit
    syncTools::syncFaceList(*this, blockedFace, orEqOp<bool>());

    regionSplit cellMark(*this, blockedFace);
    labelList cellIndi(cellMark);
    modifyRegionLabels(cellIndi);
    cRegion0Ptr_ = new labelList(cellIndi);
}


void Foam::dynamicGIBFvMesh::makeFullSnapCells() const
{
    if (fullSnapCellsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    const boolList& interfacePoints = interP();

    DynamicList<label> dlc(cells().size());

    forAll(this->cells(), cI)
    {
        const labelList& cP = this->cellPoints()[cI];
        bool foundUnsnapped = false;
        forAll(cP, pI)
        {
            const label& cPI = cP[pI];
            if (!interfacePoints[cPI])
            {
                foundUnsnapped = true;
            }
        }
        if (!foundUnsnapped)
        {
            dlc.append(cI);
        }
    }
    dlc.shrink();

    const labelList lc = labelList(dlc);

    fullSnapCellsPtr_ = new labelList(lc);
}


void Foam::dynamicGIBFvMesh::makeFullSnapCellPoints() const
{
    if (fullSnapCellPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    const labelList& fsc = fullSnapCells();
    boolList markPoints(this->points().size(), false);
    forAll(fsc, cI)
    {
        label gcI = fsc[cI];
        const labelList& cPI = this->cellPoints()[gcI];
        forAll(cPI, pI)
        {
            label gpI = cPI[pI];
            markPoints[gpI] = true;
        }
    }

    syncTools::syncPointList
    (
        *this,
        markPoints,
        plusEqOp<bool>(),
        true
    );

    DynamicList<label> dlp(points().size());
    forAll(markPoints, pI)
    {
        if (markPoints[pI])
        {
            dlp.append(pI);
        }
    }

    dlp.shrink();

    const labelList lp = labelList(dlp);

    fullSnapCellPointsPtr_ = new labelList(lp);
}


void Foam::dynamicGIBFvMesh::makePopCellPoints() const
{
    if (popCellPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popCellPointsPtr_ = new boolList(this->points().size(), false);
    boolList& popCellPoints = *popCellPointsPtr_;

    const boolList& popUpC = popUpCells();

    forAll(this->cells(), cI)
    {
        if (popUpC[cI] == 1)
        {
            const labelList& cP(this->cellPoints()[cI]);
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                popCellPoints[cPI] = true;
            }
        }
    }


    syncTools::syncPointList
    (
        *this,
        popCellPoints,
        plusEqOp<bool>(),
        true
    );

}


void Foam::dynamicGIBFvMesh::makePopSPoints() const
{
    if (popSPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popSPointsPtr_ = new boolList(this->points().size(), false);
    boolList& popSPoints = *popSPointsPtr_;

    const boolList faceIndi0 = faceIndicator0();

    const boolList& popUpC = popUpCells();

    //- mark points at the pop cells and on the interface in the previous
    // iteration
    forAll(this->cells(), cI)
    {
        if (popUpC[cI]==1)
        {
            const labelList& cP(this->cellPoints()[cI]);
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                popSPoints[cPI] = true;
            }
            const labelList& cellFaces(this->cells()[cI]);
            forAll(cellFaces, cfI)
            {
                const label& fI = cellFaces[cfI];
                //if (fI <nInternalFaces())
                {
                    if (faceIndi0[fI]==true)
                    {
                        const labelList& faceLI = this->faces()[fI];
                        forAll(faceLI, pI)
                        {
                            popSPoints[faceLI[pI]] = false;
                        }
                    }
                }
            }
        }
    }


    syncTools::syncPointList
    (
        *this,
        popSPoints,
        plusEqOp<bool>(),
        true
    );
}


void Foam::dynamicGIBFvMesh::makePopGPoints() const
{
    if (popGPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popGPointsPtr_ = new boolList(this->points().size(), false);
    boolList& popGPoints = *popGPointsPtr_;

    const boolList faceIndi = faceIndicator();

    const boolList& popUpC = popUpCells();

    forAll(this->cells(), cI)
    {
        if (popUpC[cI]==1)
        {
            const labelList& cP(this->cellPoints()[cI]);
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                popGPoints[cPI] = true;
            }
            const labelList& cellFaces(this->cells()[cI]);
            forAll(cellFaces, cfI)
            {
                const label& fI = cellFaces[cfI];
                if (faceIndi[fI]==true)
                {
                    const labelList& faceLI = this->faces()[fI];
                    forAll(faceLI, pI)
                    {
                        popGPoints[faceLI[pI]] = false;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        *this,
        popGPoints,
        plusEqOp<bool>(),
        true
    );

}


void Foam::dynamicGIBFvMesh::makeRecPoints0() const
{
    if (recPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recPoints0 = *recPoints0Ptr_;

    //const pointField& pf = this->points();

    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    const labelList& hitInd = hitIndex();
    const labelList& addrM = this->boundary()[masterId()].patch().meshPoints();

    boolList oldnewintPoints = boolList(points().size(), false);
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

    //- we want to express the point based on the coordinates of the
    //  points of the hit triangle of the stl
    //  p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a system 3x3

    forAll(addrM, pI)
    {
        label addrI = addrM[pI];
        if (!oldnewintPoints[addrI])
        {
            label hitIndexI = hitInd[addrI];
            const triFace& triList = triSM().triSurface::operator[](hitIndexI);
            //vector p = pf[addrI];
            vector p = hitPoint()[addrI];
            vector p0_1 = surP0[triList[0]];
            vector p0_2 = surP0[triList[1]];
            vector p0_3 = surP0[triList[2]];

            vector w = triList.findTriangleWeights(p, surP);

            recPoints0[addrI] = w[0]*p0_1 + w[1]*p0_2 + w[2]*p0_3;
        }
        else
        {
            recPoints0[addrI] = prevPoints_[addrI];
        }
    }
    syncPoints(recPoints0);
}


void Foam::dynamicGIBFvMesh::makeRecAllPoints0() const
{
    if (recAllPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recAllPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recAllPoints0 = *recAllPoints0Ptr_;

    pointField disp (this->points().size(), vector::zero);


//    if (!checkInterface_)
    if (true)
    {
        tmp<pointField> tsurP(triSM().points());
        const pointField& surP = tsurP();
        const pointField& surP0 = triSM().triSurface::points0();
        const labelList& hitInd = hitIndex();
        const labelList& addr = this->boundary()[masterId()].patch().meshPoints();

        if (false)
        {
            simpleVTKWriter
            (
                surP
            ).write
            (
                "stl_"+this->time().timeName()+".vtk"
            );
            simpleVTKWriter
            (
                surP0
            ).write
            (
                "stl0_"+this->time().timeName()+".vtk"
            );
        }

        //- we want to express the point based on the coordinates of the
        //  points of the hit triangle of the stl
        //  p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a system 3x3


        //- first for the master points

        boolList pointMoved(this->points().size(), false);

        forAll(addr, pI)
        {
            label addrI = addr[pI];
            label hitIndexI = hitInd[addrI];
            const triFace& triList = triSM().triSurface::operator[](hitIndexI);
            vector p = hitPoint()[addrI];
            vector p0_0 = surP0[triList[0]];
            vector p0_1 = surP0[triList[1]];
            vector p0_2 = surP0[triList[2]];

/*
            vector p_0 = surP[triList[0]];
            vector p_1 = surP[triList[1]];
            vector p_2 = surP[triList[2]];
*/

            vector w = triList.findTriangleWeights(p, surP);
            vector newLoc = w[0]*p0_0 + w[1]*p0_1 + w[2]*p0_2;
            disp[addrI] = newLoc - recAllPoints0[addrI];
            pointMoved[addrI] = true;
        }

        syncTools::syncPointList
        (
            *this,
            pointMoved,
            plusEqOp<bool>(),
            true
        );


        const labelList& addrS = this->boundary()[slaveId()].patch().meshPoints();

        //- we want to express the point based on the coordinates of the
        //  points of the hit triangle of the stl
        //  p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a system 3x3


        //- second the points in slave patch that they have not been updated
        //- from the master loop points

        forAll(addrS, pI)
        {
            label addrI = addrS[pI];
            if (pointMoved[addrI] == false)
            {
                label hitIndexI = hitInd[addrI];
                const triFace& triList = triSM().triSurface::operator[](hitIndexI);
                vector p = hitPoint()[addrI];
                vector p0_1 = surP0[triList[0]];
                vector p0_2 = surP0[triList[1]];
                vector p0_3 = surP0[triList[2]];

                vector w = triList.findTriangleWeights(p, surP);

                vector newLoc = w[0]*p0_1 + w[1]*p0_2 + w[2]*p0_3;
                disp[addrI] = newLoc  - recAllPoints0[addrI];
                pointMoved[addrI] = true;
            }
        }

        const labelList& fscp = fullSnapCellPoints();

        syncTools::syncPointList
        (
            *this,
            disp,
            maxMagEqOp(),
            vector::zero
        );
        recAllPoints0 += disp;

        const boolList& interPoints = interP();
        const boolList& interPoints0 = interP0();

        computeOldPositionsInUnsnappedCells(recAllPoints0, fscp);

        const boolList& meshQualityP = meshQualityPoints();

        forAll(meshQualityP, pI)
        {
            if
            (
                meshQualityP[pI] && interPoints[pI] && interPoints0[pI]
            )
            {
                recAllPoints0[pI]  = prevPoints_[pI];
            }
        }

        /*
        DynamicList<label> dFaces(faces().size());
        if (Pstream::myProcNo()==1)
        {
            label fNo = 265589;

            forAll(faces()[fNo], pI)
            {
                label gpI = faces()[fNo][pI];
                Pout<< gpI <<  tab << recAllPoints0[gpI] << tab << points()[gpI] << endl;
            }
            dFaces.append(fNo);
        }
        dFaces.shrink();
        simpleVTKWriter
        (
            this->faces(),
            dFaces,
            recAllPointsBU
        ).write
        (
            "face0BU_"+this->time().timeName()+".vtk"
        );
        simpleVTKWriter
        (
            this->faces(),
            dFaces,
            recAllPoints0
        ).write
        (
            "face0_"+this->time().timeName()+".vtk"
        );
        simpleVTKWriter
        (
            this->faces(),
            dFaces,
            points()
        ).write
        (
            "face_"+this->time().timeName()+".vtk"
        );
        simpleVTKWriter
        (
            this->faces(),
            fl(),
            points()
        ).write
        (
            "inject"+this->time().timeName()+".vtk"
        );
    */

    }
    else
    {
        recAllPoints0 = prevPoints_;
    }
}


void Foam::dynamicGIBFvMesh::makePopUpCells() const
{
    const labelList& cReg = cRegion();
    const labelList& cReg0 = cRegion0();

    popUpCellsPtr_ = new boolList(this->cells().size(), false);
    boolList& popUpCells = *popUpCellsPtr_;

    forAll(this->cells(), cI)
    {
        if (cReg[cI] != cReg0[cI])
        {
            popUpCells[cI] = 1;
        }
        else
        {
            popUpCells[cI] = 0;
        }
    }
    if (time().outputTime() && debugMode_)
    {
        writeScalarField("popCells", popUpCells);
    }
}



void Foam::dynamicGIBFvMesh::makeBoundaryPoints() const
{
    if (boundaryPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boundaryPointsPtr_ = new boolList(this->points().size(), false);

    boolList& boundaryPoints = *boundaryPointsPtr_;

    for (label fI = this->nInternalFaces(); fI < this->nFaces(); fI++)
    {
        label patchI = this->boundaryMesh().whichPatch(fI);
        if
        (
            !isA<indirectPolyPatch>(this->boundaryMesh()[patchI]) &&
            !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
            !isA<emptyFvPatch>(this->boundary()[patchI]) &&
            !this->boundary()[patchI].coupled()
        )
        {
            const labelList& faceI = this->faces()[fI];
            forAll(faceI, fII)
            {
                boundaryPoints[faceI[fII]] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        *this,
        boundaryPoints,
        plusEqOp<bool>(),
        true
    );
}


void Foam::dynamicGIBFvMesh::extendBoundaryBaseMesh
(
    const labelList& interestFaces,
    vectorField& neiCC
) const
{
    const vectorField& baseCf = *baseCf_;
    const vectorField& baseCC = *baseCC_;
    forAll(interestFaces, fI)
    {
        label interestFacesI = interestFaces[fI];
        if (interestFacesI >= this->nInternalFaces())
        {
            label pI = this->boundaryMesh().whichPatch(interestFacesI);
            const polyPatch& pp = this->boundary()[pI].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                if (!pp.coupled())
                {
                    const label& on = this->faceOwner()[interestFacesI];
                    const vector& ccI = baseCC[on];

                    vector vfc = baseCf[interestFacesI]-ccI;
                    scalar mvfc = mag(vfc);
                    const vector& nf = (vfc)/(mvfc+SMALL);

                    //- extend base mesh face center points by cfCC vector
                    neiCC[fI] = baseCf[interestFacesI]+(nf*mvfc);
                }
            }
        }
    }
    if (false)
    {
        simpleVTKWriter
        (
            neiCC
        ).write
        (
            "extBoundary_"+this->time().timeName()+".vtk"
        );
    }
}


faceList Foam::dynamicGIBFvMesh::preparePatch
(
    const fvPatch& gibPatch
)
{
    faceList faces(gibPatch.patch().localFaces());

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    const boolList& flipmap = gibPolyPatch.fm();

    forAll(faces, fI)
    {
        if (flipmap[fI])
        {
            faces[fI].flip();
        }
    }

    return faces;
}


void Foam::dynamicGIBFvMesh::syncPoints
(
    pointField& points
) const
{
    //--- syncing the processor boundary points ---//
    pointField np1 (points);
    pointField np2 (points);

    syncTools::syncPointPositions
    (
        *this,
        np1,
        minEqOp<point>(),           // combine op
        point(GREAT,GREAT,GREAT)    // null
    );
    syncTools::syncPointPositions
    (
        *this,
        np2,
        maxEqOp<point>(),           // combine op
        point(-GREAT,-GREAT,-GREAT)    // null
    );

    forAll(points, pI)
    {
        points[pI] = (np1[pI]+np2[pI])/2;
    }
}


labelList Foam::dynamicGIBFvMesh::checkingIntersectingFaces() const
{
    const labelList& fl = fl0();

    boolList interestedFace(this->nFaces(), false);
    boolList constraintFace(this->nFaces(), false);

    //- if initiallization take all the faces
    label ggibSize = fl.size();

    reduce(ggibSize, sumOp<label>());

    // mark all the faces we dont want to be included
    forAll(this->boundary(), pI)
    {
        const polyPatch& pp = this->boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if
            (
                isA<wedgeFvPatch>(this->boundary()[pI]) ||
                isA<emptyFvPatch>(this->boundary()[pI])
            )
            {
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    constraintFace[gfI] = true;
                }
            }
        }
    }

    if (ggibSize == 0)
    {
        //-if old time doens't exist (for example when runing createGIB)
        // add all the faces except the constraint faces
        forAll(interestedFace, fI)
        {
            if (constraintFace[fI] == false)
            {
                interestedFace[fI] = true;
            }
        }
        forAll(this->boundary(), pI)
        {
            const polyPatch& pp = this->boundary()[pI].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                if ((!pp.coupled()) && !includeWalls())
                {
                    forAll(pp, pfI)
                    {
                        label gfI = pp.start() + pfI;
                        interestedFace[gfI] = false;
                    }
                }
            }
        }
        DynamicList<label> dfl(this->faces().size());

        forAll(interestedFace, fI)
        {
            if (interestedFace[fI])
            {
                dfl.append(fI);
            }
        }
        dfl.shrink();
        return labelList(dfl);
    }


    // mark all the faces of the previous interface
    forAll(fl, fI)
    {
        label flI = fl[fI];
        interestedFace[flI] = true;
    }

    //- mark all the boundaries
    forAll(this->boundary(), pI)
    {
        const polyPatch& pp = this->boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            forAll(pp, pfI)
            {
                label gfI = pp.start() + pfI;
                if (!constraintFace[gfI])
                {
                    interestedFace[gfI] = true;
                }
            }
        }
    }


    boolList bandFaces(this->nFaces(), false);

    for (int i = 1; i<4; i++) //controling the band
    {
        //- mark all the face-cell-faces of the internal faces
        forAll(this->neighbour(), fI)
        {
            if (interestedFace[fI])
            {
                const label& on = this->owner()[fI];
                const label& nb = this->neighbour()[fI];
                {
                    const labelList& cellFaces = this->cells()[on];
                    forAll(cellFaces, cfI)
                    {
                        label cfII = cellFaces[cfI];
                        if (!constraintFace[cfII])
                        {
                            bandFaces[cfII] = true;
                        }
                    }
                }
                {
                    const labelList& cellFaces = this->cells()[nb];
                    forAll(cellFaces, cfI)
                    {
                        label cfII = cellFaces[cfI];
                        if (!constraintFace[cfII])
                        {
                            bandFaces[cfII] = true;
                        }
                    }
                }
            }
        }

        //- mark all the face-cell-faces of the boundaryCells
        forAll(this->boundary(), pI)
        {
            const polyPatch& pp = this->boundary()[pI].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                const labelList& fcs = this->boundary()[pI].faceCells();
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    if (interestedFace[gfI])
                    {
                        label fc = fcs[pfI];
                        const labelList& cellFaces = this->cells()[fc];
                        forAll(cellFaces, cfI)
                        {
                            label cfII = cellFaces[cfI];
                            if (!constraintFace[cfII])
                            {
                                bandFaces[cfII] = true;
                            }
                        }
                    }
                }
            }
        }

        forAll(interestedFace, fI)
        {
            if (bandFaces[fI])
            {
                interestedFace[fI] = true;
            }
        }
    }

    forAll(this->boundary(), pI)
    {
        const polyPatch& pp = this->boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if ((!pp.coupled()) && !includeWalls())
            {
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    interestedFace[gfI] = false;
                }
            }
        }
    }

    DynamicList<label> dfl(this->faces().size());

    forAll(interestedFace, fI)
    {
        if (interestedFace[fI])
        {
            dfl.append(fI);
        }
    }
    dfl.shrink();

    return labelList(dfl);
}


bool Foam::dynamicGIBFvMesh::checkIfEdgeToleranceIntersecting
(
    const face& face1,
    const face& face2
) const
{
    int samePoints = 0;
    bool tolIssue = false;

    forAll(face1, f1)
    {
        forAll(face2, f2)
        {
            if (face1[f1] == face2[f2])
            {
                samePoints +=1;
            }
        }
    }
    if (samePoints==2)
    {
        tolIssue = true;
    }

    return tolIssue;
}


void Foam::dynamicGIBFvMesh::doTangentialBoundaryMotion
(
    pointField& newPoints
) const
{
    const vectorField& baseSf = *baseSf_;
    const pointField& basePoints = *basePoints_;

    const boolList& bPoints = boundaryPoints();

    const labelList& fList = fl();

    boolList gibPoints = boolList(basePoints.size(), false);

    forAll(fList, fI)
    {
        label fII = fList[fI];
        const face& faceI = this->faces()[fII];
        forAll(faceI, pI)
        {
            gibPoints[faceI[pI]] = true;
        }
    }

    boolList gibBoundaryPoints = boolList(basePoints.size(), false);

    forAll(gibBoundaryPoints, pI)
    {
        if (bPoints[pI] && gibPoints[pI])
        {
            gibBoundaryPoints[pI] = true;
        }
    }

    const labelListList& pFaces = this->pointFaces();
    forAll(gibBoundaryPoints, pI)
    {
        if (gibBoundaryPoints[pI])
        {
            vector dis = newPoints[pI] - basePoints[pI];

            const labelList& pFaceI = pFaces[pI];
            forAll(pFaceI, fI)
            {
                label faceI = pFaceI[fI];
                if (faceI >= this->nInternalFaces())
                {
                    label patchI = this->boundaryMesh().whichPatch(faceI);
                    const polyPatch& poly = this->boundaryMesh()[patchI];
                    if (isA<directPolyPatch>(poly))
                    {
                        if
                        (
                            !(
                                this->boundaryMesh()[patchI].coupled() ||
                                isA<wedgePolyPatch>(this->boundaryMesh()[patchI]) ||
                                isA<emptyPolyPatch>(this->boundaryMesh()[patchI])
                            )
                        )
                        {
                            vector n = baseSf[faceI]/mag(baseSf[faceI]);
                            dis -= (dis&n)*n;
                        }
                    }
                }
            }
            newPoints[pI] = basePoints[pI] + dis;
        }
    }

}




tmp<Foam::pointField> Foam::dynamicGIBFvMesh::findSnappedPoints
(
    const bool& fromBase
) const
{
    const pointField& baseP(*basePoints_);

    //- Initialize the snapped with the current points
    tmp<pointField> npt(new pointField(this->points()));
    pointField& np = npt.ref();
    if (fromBase)
    {
        np = baseP;
    }
    else
    {
        np = this->points();
    }

    const boolList& interPoints = interP();
    const boolList& interPoints0 = interP0();
    DynamicList<label> dIntPoints(this->points().size());

    forAll(interPoints, pI)
    {
        if (interPoints[pI])
        {
            dIntPoints.append(pI);
        }
    }
    dIntPoints.shrink();

    if (!fromBase)
    {
        //  The points of the current GIB are reverted to the basePoint location
        forAll(dIntPoints, pI)
        {
            const label gpI = dIntPoints[pI];
            np[gpI] = baseP[gpI];
        }
    }
    pointField dIntPointLoc = pointField(dIntPoints.size(), vector::zero);

    List<pointIndexHit> nearest;

    boolList multInterPoints = findMultInterPoints();

    const boolList& popCellP = popCellPoints();

    //- mark points which belong to the cells that all their points are snapped
    const labelList& fscp = fullSnapCellPoints();
    boolList gfullSnapCellPoints(this->points().size(), false);
    forAll(fscp, pI)
    {
        label gpI = fscp[pI];
        gfullSnapCellPoints[gpI] = true;
    }

    const boolList& bPoints = boundaryPoints();

    forAll(dIntPoints, pI)
    {
        label gpI = dIntPoints[pI];
        dIntPointLoc[pI] = baseP[gpI];
    }

    smoothInternalBasePoints(dIntPointLoc, dIntPoints);

    forAll(dIntPoints, pI)
    {
        label gpI = dIntPoints[pI];
        if
        (
            (
                (
                    !interPoints0[gpI]
                    ||  multInterPoints[gpI]
                    ||  gfullSnapCellPoints[gpI]
                    ||  popCellP[gpI]
                )
                && !bPoints[gpI]
            )
            || fromBase
        )
        {
        }
        else
        {
            dIntPointLoc[pI] = this->points()[gpI];
        }
    }


    triSM().findNearest
    (
        dIntPointLoc,
        scalarField(dIntPointLoc.size(), sqr(GREAT)),
        nearest
    );

    labelIOList& hitIndex = *hitIndexPtr_;
    hitIndex.resize(this->points().size(), -1);

    hitPointPtr_ = new vectorField(this->points().size(), vector::zero);
    vectorField& hitPoint = *hitPointPtr_;

    forAll(dIntPoints, pI)
    {
        const pointIndexHit& nearI = nearest[pI];
        np[dIntPoints[pI]] = nearI.hitPoint();
        hitIndex[dIntPoints[pI]] = nearI.index();
    }

    syncPoints(np);
    hitPoint = np;

    pointField fromPoints(this->points().size(), vector::zero);
    forAll(dIntPoints, pI)
    {
        label gpI = dIntPoints[pI];
        fromPoints[gpI] = dIntPointLoc[pI];
    }


    //checkConcaveBoundaryPoints(baseP, np);
    fullSnappedPointsTreatment(np, baseP);

    unsnapForQuality(np, baseP);

    syncPoints(np);

    return npt;
}


void Foam::dynamicGIBFvMesh::unsnapForQuality
(
    pointField& newPoints,
    const pointField& fromPoints
) const
{
    pointField snappedPoints = newPoints;
    const boolList& isInterfacePoint = interP();
    scalarField r(newPoints.size(), 0);
    const label maxUnsnapIntervals = 10;
    label unsnapStep;

    if (meshQualityPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    meshQualityPointsPtr_ = new boolList(points().size(), false);
    boolList& meshQualityPoints = *meshQualityPointsPtr_;

    boolList hasInterPoint(this->faces().size(), false);
    forAll(hasInterPoint, fi)
    {
        const face& f = this->faces()[fi];
        forAll(f, fpi)
        {
            if (isInterfacePoint[f[fpi]])
            {
                hasInterPoint[fi] = true;
                break;
            }
        }
    }

    for (unsnapStep = 0; unsnapStep < maxUnsnapIntervals; unsnapStep++)
    {
        bool qualityOK = true;
        scalarList nonOrth(this->faces().size(), GREAT);
        scalarList ownPyrVol(this->faces().size(), scalar(0));
        scalarList neiPyrVol(this->faces().size(), scalar(0));

        vectorField d(this->faces().size(), vector::zero);
        forAll(this->faces(), fi)
        {
            if (hasInterPoint[fi])
            {
                label celli = this->faceOwner()[fi];
                const face& f = this->faces()[fi];
                point ownCentre
                (
                    this->cells()[celli].centre(newPoints, this->faces())
                );

                // Calc cell-centre-to-cell-centre vectors
                d[fi] = -ownCentre;

                // Calculate pyramid volumes
                point fc(f.centre(newPoints));
                vector fa(f.areaNormal(newPoints));
                ownPyrVol[fi] =
                   -primitiveMeshTools::pyramidVol(ownCentre, fc, fa);

                if (fi < this->faceNeighbour().size())
                {
                    label cellj = this->faceNeighbour()[fi];
                    point neiCentre
                    (
                        this->cells()[cellj].centre(newPoints, this->faces())
                    );

                    d[fi] += neiCentre;

                    neiPyrVol[fi] =
                        primitiveMeshTools::pyramidVol(neiCentre, fc, fa);
                }
                else
                {
                    // This gets cancelled for coupled boundaries
                    d[fi] += f.centre(newPoints);
                }
            }
        }

        // Calculate across coupled boundaries
        syncTools::syncFaceList(*this, d, minusEqOp<vector>());

        forAll(this->faces(), fi)
        {
            if (hasInterPoint[fi])
            {
                // Calc non-orth
                const face& f = this->faces()[fi];
                vector nf =
                    f.areaNormal(newPoints)/stabilise(f.mag(newPoints), SMALL);
                nonOrth[fi] = (nf&d[fi])/stabilise(mag(d[fi]), SMALL);
            }
        }

        boolList unsnapP(newPoints.size(), false);
        forAll(this->faces(), fi)
        {
            if (hasInterPoint[fi])
            {
                if
                (
                    nonOrth[fi] < maxNonOrth_
                 || ownPyrVol[fi] < 0
                 || neiPyrVol[fi] < 0
                )
                {
                    qualityOK = false;
                    const face& f = this->faces()[fi];
                    forAll(f, fpi)
                    {
                        unsnapP[f[fpi]] = true;
                    }
                }
            }
        }

        reduce(qualityOK, andOp<bool>());
        if (qualityOK)
        {
            break;
        }

        // Sync here as unsnapP might have got changed by a non-coupled face
        // on one side
        syncTools::syncPointList(*this, unsnapP, orEqOp<bool>(), false);
        forAll(meshQualityPoints, pi)
        {
            meshQualityPoints[pi] = meshQualityPoints[pi] || unsnapP[pi];
        }

        forAll(unsnapP, pi)
        {
            if (unsnapP[pi])
            {
                r[pi] = scalar(unsnapStep+1)/(maxUnsnapIntervals);
                newPoints[pi] =
                    (
                        snappedPoints[pi]
                      - (snappedPoints[pi]-fromPoints[pi])*r[pi]
                    );
            }
        }
    }
    if (debugMode_)
    {
        scalar totalr(0), maxr(0);
        label numIntP(0);
        forAll(r, pi)
        {
            if (isInterfacePoint[pi])
            {
                totalr += r[pi];
                maxr = max(r[pi], maxr);
                numIntP++;
            }
        }
        Foam::reduce(
            std::tie(totalr, maxr, numIntP),
            ParallelOp<sumOp<scalar>, maxOp<scalar>, sumOp<label>>{},
            comm()
        );

        if (!numIntP)
        {
            // Avoid /0 below
            numIntP = 1;
        }
        Info<< "GIB mesh-quality constraint: Unsnapping fraction  "
            << "Mean = " << totalr/numIntP << ", max = " << maxr << endl;
    }
}


void Foam::dynamicGIBFvMesh::fullSnappedPointsTreatment
(
    pointField& newPoints,
    const pointField& fromPoints
) const
{
    const labelList& fscp = fullSnapCellPoints();

    forAll(fscp, pI)
    {
        label gpI = fscp[pI];
        newPoints[gpI] -= unsnapVar_*(newPoints[gpI] - fromPoints[gpI]);
    }
}

void Foam::dynamicGIBFvMesh::computeOldPositionsInUnsnappedCells
(
    pointField& recAllPoints0,
    const labelList& fscp
) const
{
    forAll(fscp, pI)
    {
        label gpI = fscp[pI];
        recAllPoints0[gpI]  = this->points()[gpI];
    }
}


label Foam::dynamicGIBFvMesh::gibFaceZone() const
{
    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>
        (
            this->boundary()[masterId()].patch()
        );

    return gibPolyPatch.zoneId();
}


label Foam::dynamicGIBFvMesh::findOrCreateCellZone() const
{
    const label& fZoneId = gibFaceZone();
    const word faceZoneName = this->faceZones()[fZoneId].name();

    const word cellZoneName = "inactive_"+faceZoneName;
    const label& cZoneId = this->cellZones().findZoneID(cellZoneName);

    if (cZoneId != -1)
    {
        return cZoneId;
    }
    else
    {
        Info<< "cellZone was not found:" << endl;
        Info<< "Constructing new cellZone called " << cellZoneName << endl;
        const polyMesh& pMesh = *this;
        polyMesh& cpMesh = const_cast<polyMesh&>(pMesh);
        label oldZoneSize = cpMesh.cellZones().size();
        cpMesh.cellZones().setSize(oldZoneSize + 1);

        label index = oldZoneSize;
        cpMesh.cellZones().set
        (
            index,
            new cellZone
            (
                cellZoneName,
                labelUList(),
                index,
                cpMesh.cellZones()
            )
        );

        return index;
    }
}

void Foam::dynamicGIBFvMesh::calculateFlipmap() const
{
    updateSolidCellZone();
    const cellZone& cZone = this->cellZones()[findOrCreateCellZone()];
    boolList markCells = boolList(this->cells().size(), false);
    const labelList& fList = fl();
    boolList fm = boolList(fList.size(), false);

    forAll(cZone, cI)
    {
        const label& czI = cZone[cI];
        markCells[czI] = true;
    }

    forAll(fList, fI)
    {
        const label& fII = fList[fI];
        {
            const label& on = this->faceOwner()[fII];
            if (markCells[on])
            {
                fm[fI] = true;
            }
        }
    }
    fmPtr_ = new boolList (fm);
}

void Foam::dynamicGIBFvMesh::regionVisDebug() const
{
    const labelList& cReg = cRegion();
    volScalarField cellRegionF
    (
        IOobject
        (
            "cellRegionF",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("0", dimless, 0),
        "calculated"
    );
    forAll(cellRegionF, cI)
    {
        cellRegionF[cI] = cReg[cI];
    }
    if (time().outputTime())
    {
        cellRegionF.write();
    }
}


boolList Foam::dynamicGIBFvMesh::findMultInterPoints() const
{
    boolList multInterPoints(this->points().size(), false);
    const boolList& multInterFaces = *multInterFacesPtr_;
    forAll(multInterFaces, fI)
    {
        if (multInterFaces[fI] == true)
        {
            const labelList& fpI = this->faces()[fI];
            forAll(fpI, pI)
            {
                multInterPoints[fpI[pI]] = true;
            }
        }
    }

    return multInterPoints;
}



void Foam::dynamicGIBFvMesh::popShrinkFields()
{

    const boolList& popSP = popSPoints();

    //- mark all the popSP neightbour faces & points
    // iteration

    boolList fPopIndi  = boolList(this->faces().size(), false);
    boolList popIndi  = boolList(this->cells().size(), false);
    forAll(this->points(), pI)
    {
        const labelList& pFacesI = this->pointFaces()[pI];
        const labelList& pCellsI = this->pointCells()[pI];
        if (popSP[pI])
        {
            forAll(pFacesI, fI)
            {
                const label& pFI = pFacesI[fI];
                if (pFI<this->nInternalFaces())
                {
                    fPopIndi[pFI] = true;
                }
                else
                {
                    label patchI = this->boundaryMesh().whichPatch(pFI);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                        !isA<emptyFvPatch>(this->boundary()[patchI])
                    )
                    {
                        fPopIndi[pFI] = true;
                    }
                }
            }
            forAll(pCellsI, cI)
            {
                const label& pCI = pCellsI[cI];
                popIndi[pCI] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, fPopIndi, orEqOp<bool>());

    DynamicList<label> fPopD(this->faces().size());
    forAll(this->faces(), fI)
    {
        //if(fI<this->nInternalFaces())
        {
            if (fPopIndi[fI])
            {
                fPopD.append(fI);
            }
        }
    }
    fPopD.shrink();
    labelList fPop(fPopD);

    const surfaceScalarField& phiPop = phiTrS();

    PopFields<scalar, fvPatchField, volMesh>(fPop, phiPop, popIndi, false);
    PopFields<vector, fvPatchField, volMesh>(fPop, phiPop, popIndi, false);
    PopFields<sphericalTensor, fvPatchField, volMesh>(fPop, phiPop, popIndi, false);
    PopFields<symmTensor, fvPatchField, volMesh>(fPop, phiPop, popIndi, false);
    PopFields<tensor, fvPatchField, volMesh>(fPop, phiPop, popIndi, false);
}


void Foam::dynamicGIBFvMesh::popGrowFields()
{
    const scalarField& cV0 = this->V0();

    const boolList& popGP = popGPoints();

    const boolList& popUpC = popUpCells();

    boolList fPopIndi  = boolList(this->faces().size(), false);
    boolList popIndi  = boolList(this->cells().size(), false);
    forAll(this->points(), pI)
    {
        const labelList& pFacesI = this->pointFaces()[pI];
        const labelList& pCellsI = this->pointCells()[pI];
        if (popGP[pI])
        {
            forAll(pFacesI, fI)
            {
                const label& pFI = pFacesI[fI];
                if (pFI<this->nInternalFaces())
                {
                    fPopIndi[pFI] = true;
                }
                else
                {
                    label patchI = this->boundaryMesh().whichPatch(pFI);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                        !isA<emptyFvPatch>(this->boundary()[patchI])
                    )
                    {
                        fPopIndi[pFI] = true;
                    }
                }
            }
            forAll(pCellsI, cI)
            {
                const label& pCI = pCellsI[cI];
                popIndi[pCI] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, fPopIndi, orEqOp<bool>());

    DynamicList<label> fPopD(this->faces().size());
    forAll(this->faces(), fI)
    {
        if (fPopIndi[fI])
        {
            fPopD.append(fI);
        }
    }
    fPopD.shrink();
    labelList fPop(fPopD);

    const surfaceScalarField& phiPop = phiTrG();

    scalarField sV = scalarField(this->cells().size(), 0.0);

    //- oldVgrow -> volume when we start growing from 0 sized volume
    scalarField oldVolgrow = cV0;


    {
        forAll(fPop, fI)
        {
            const label& fPopI = fPop[fI];
            if (fPopI<this->nInternalFaces())
            {
                scalar sVf = phiPop[fPopI]*this->time().deltaTValue();

                const label& on = this->owner()[fPopI];
                const label& nb = this->neighbour()[fPopI];
                sV[on] += sVf;
                sV[nb] -= sVf;
            }
            else
            {
                const label& on = this->faceOwner()[fPopI];
                label patchI = this->boundaryMesh().whichPatch(fPopI);
                label lpfI = fPopI - this->boundaryMesh()[patchI].start();
                scalar sVf = phiPop.boundaryField()[patchI][lpfI]*
                    this->time().deltaTValue();
                sV[on] += sVf;
            }
        }

        forAll(popIndi, cI)
        {
            if (popIndi[cI])
            {
                oldVolgrow[cI] -= sV[cI];
            }
        }
    }

    forAll(this->cells(), cI)
    {
        if (popUpC[cI]==1)
        {
            oldVolgrow[cI] = 0;
        }
    }

    PopFields2<scalar, fvPatchField, volMesh>(fPop, phiPop, popIndi, oldVolgrow);
    PopFields2<vector, fvPatchField, volMesh>(fPop, phiPop, popIndi, oldVolgrow);
    PopFields2<sphericalTensor, fvPatchField, volMesh>
    (
        fPop, phiPop, popIndi, oldVolgrow
    );
    PopFields2<symmTensor, fvPatchField, volMesh>(fPop, phiPop, popIndi, oldVolgrow);
    PopFields2<tensor, fvPatchField, volMesh>(fPop, phiPop, popIndi, oldVolgrow);
}


void Foam::dynamicGIBFvMesh::correctV0()
{
    scalarField& V =
        const_cast<scalarField&>(dynamic_cast<const scalarField&>(this->V()));
    scalarField& V0 =
        const_cast<scalarField&>(dynamic_cast<const scalarField&>(this->V0()));
    surfaceScalarField& phiB = const_cast<surfaceScalarField&>(this->phiB());
    dimensionedScalar deltaT = this->time().deltaT();

    // Recompute V0
    V0 = V*(1 - scalarField(fvc::div(phiB))*deltaT.value());

    // Make negative V0 to small positive
    bool V0Negative = false;
    forAll(V0, celli)
    {
        if (V0[celli] < SMALL)
        {
            V0[celli] = SMALL;
            V0Negative = true;
        }
    }

    // The flag has to be synchronised otherwise in the equation solution for
    // correction the boundaries can't be synced and solver will freeze
    reduce(V0Negative, orOp<bool>());

    // Correct GIB fluxes to correspond to non-negative V0
    // This happens for the badly snapped cells. Sometimes it can be avoided
    // by reducing the timestep.
    if (V0Negative)
    {
        // This function replaces modifyVolSqueezedCells
        // and fixes the issue with volumes not coressponding to fluxes.
        // However, that does require solving equation (adds on computational time)
        negativeV0CorrectGIBFluxes(phiB, V, V0, deltaT);
    }
}


void Foam::dynamicGIBFvMesh::negativeV0CorrectGIBFluxes
(
    surfaceScalarField& phiB,
    const scalarField& V,
    const scalarField& V0,
    const dimensionedScalar& deltaT
)
{
    // Initialize BCs list for GIBMeshPhiCorr to zero-gradient
    wordList GIBMeshPhiCorrTypes
    (
        phiB.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of GIBMeshPhiCorr to fixed-value for patches at which p is fixed
    forAll(phiB.boundaryField(), patchi)
    {
        if (phiB.boundaryField()[patchi].fixesValue())
        {
            GIBMeshPhiCorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField GIBMeshPhiCorr
    (
        IOobject
        (
            "GIBMeshPhiCorr",
            this->time().timeName(),
            *this
        ),
        *this,
        dimensionedScalar
        (
            "GIBMeshPhiCorr",
            phiB.dimensions()*dimTime/dimLength,
            0.0
        ),
        GIBMeshPhiCorrTypes
    );
    this->schemes().setFluxRequired(GIBMeshPhiCorr.name());

    volScalarField divVolRequired
    (
        IOobject
        (
            "divVolRequired",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimensionedScalar("divVolRequired", dimless/dimTime, 0.0)
    );
    divVolRequired.primitiveFieldRef() = ((V - V0)/(V*deltaT.value()));

    fvScalarMatrix phiCorrEqn
    (
        fvm::laplacian(1.0/deltaT, GIBMeshPhiCorr)
     ==
        fvc::div(phiB) - divVolRequired
    );

    // Solver settings from cellDisplacement fvSolution sub-dictionary
    phiCorrEqn.solve(this->solution().solver("cellDisplacement"));
    phiB -= phiCorrEqn.flux();
}


void Foam::dynamicGIBFvMesh::modifyVolSqueezedCells() const
{
    const scalarField& cV0 = this->V0();
    const scalarField& cV = this->V();

    DynamicList<label> dcV0(cells().size());
    DynamicList<label> dcV(cells().size());

    forAll(cells(), cI)
    {
        if (cV0[cI] < SMALL)
        {
            dcV0.append(cI);
        }
        if (cV[cI] < SMALL)
        {
            dcV.append(cI);
        }
    }
    dcV0.shrink();
    dcV.shrink();


    scalarField& V0 =  const_cast<scalarField&>(cV0);
    scalarField& V =  const_cast<scalarField&>(cV);

    const surfaceScalarField& cmPhi = this->phi();
    surfaceScalarField& mPhi =  const_cast<surfaceScalarField&>(cmPhi);
    surfaceScalarField::Boundary& mPhibf = mPhi.boundaryFieldRef();

    forAll(dcV0, cI)
    {
        label gcI = dcV0[cI];

        if (V[gcI] < SMALL)
        {
            V[gcI] = SMALL;
        }
        V0[gcI] = V[gcI];

        const labelList& cellFaces(this->cells()[gcI]);
        forAll(cellFaces, cfI)
        {
            label faceI = cellFaces[cfI];
            if (faceI < nInternalFaces())
            {
                mPhi[faceI] = 0;
            }
            else
            {
                label pI = this->boundaryMesh().whichPatch(faceI);
                const polyPatch& pp = this->boundary()[pI].patch();
                if (mPhibf[pI].size() != 0)
                {
                    if (!isA<indirectPolyPatch>(pp))
                    {
                        label pfI = faceI - pp.start();
                        mPhibf[pI][pfI] = 0;
                    }
                }
            }
        }
    }

    forAll(V, gcI)
    {
        if (V[gcI] < SMALL)
        {
            V[gcI] = SMALL;
            V0[gcI] = SMALL;
            const labelList& cellFaces(this->cells()[gcI]);
            forAll(cellFaces, cfI)
            {
                label faceI = cellFaces[cfI];
                if (faceI < nInternalFaces())
                {
                    mPhi[faceI] = 0;
                }
                else
                {
                    label pI = this->boundaryMesh().whichPatch(faceI);
                    const polyPatch& pp = this->boundary()[pI].patch();
                    if (mPhibf[pI].size() != 0)
                    {
                        if (!isA<indirectPolyPatch>(pp))
                        {
                            label pfI = faceI - pp.start();
                            mPhibf[pI][pfI] = 0;
                        }
                    }
                }
            }
        }
    }

    if (debugMode_)
    {

        if (dcV0.size())
        {
            OFstream strV0
            (
                this->time().path()/
                "cells_V0_"+this->time().timeName()+".obj"
            );

            meshTools::writeOBJ
            (
                strV0,
                this->cells(),
                this->faces(),
                this->points(),
                labelList(dcV0)
            );
        }

        if (dcV.size())
        {
            OFstream strV
            (
                this->time().path()/
                "cells_V_"+this->time().timeName()+".obj"
            );
            meshTools::writeOBJ
            (
                strV,
                this->cells(),
                this->faces(),
                this->points(),
                labelList(dcV)
            );
        }
    }
}


void Foam::dynamicGIBFvMesh::resetMeshFluxes()
{
    //- used to cancel the meshPhi with phi at boundary
    //- to make the relative 0, but making that the mesh fluxes are inconsistent
    const surfaceScalarField& phiBoundary = phiB();
    const surfaceScalarField& cmPhi = this->phi();
    surfaceScalarField& mPhi =  const_cast<surfaceScalarField&>(cmPhi);
    mPhi = phiBoundary;
    volVectorField& U = lookupObjectRef<volVectorField>("U");

    U.correctBoundaryConditions();
    U.oldTime().correctBoundaryConditions(); //- not consistent
}

void Foam::dynamicGIBFvMesh::correctVelocityFlux()
{
    surfaceScalarField& phi =
        const_cast<surfaceScalarField&>(lookupObject<surfaceScalarField>("phi"));
    const volVectorField& U = this->lookupObject<volVectorField>("U");

    phi = fvc::interpolate(U) & this->Sf();
}


void Foam::dynamicGIBFvMesh::visCells()
{
    masterFCells_ = 0;
    slaveFCells_ = 0;

    const labelList& fc1 = this->boundary()[masterGIB_].faceCells();
    const labelList& fc2 = this->boundary()[slaveGIB_].faceCells();
    forAll(fc1, fcI)
    {
        const label& fc1I = fc1[fcI];
        if (fc1I<this->nInternalFaces())
        {
            masterFCells_[fc1I] = 1;
        }
        else
        {
            masterFCells_[fc1I] = 1;
        }
    }
    forAll(fc2, fcI)
    {
        const label& fc2I = fc2[fcI];
        if (fc2I<this->nInternalFaces())
        {
            slaveFCells_[fc2I] = 1;
        }
        else
        {
            slaveFCells_[fc2I] = 1;
        }
    }
}


void Foam::dynamicGIBFvMesh::faceCellsVisDebug()
{
    if (debugMode_)
    {
        visCells();
        writeScalarField("masterCells", masterFCells_);
        writeScalarField("slaveCells", slaveFCells_);
    }

    writeScalarField("cRegion", cRegion());
}


void Foam::dynamicGIBFvMesh::findGIBPatches()
{
    patchGIBPtr_ = new boolList(this->boundary().size(), false);
    boolList& patchGIB = *patchGIBPtr_;
    forAll(this->boundary(), pI)
    {
        const polyPatch& poly = this->boundary()[pI].patch();
        if (isA<indirectPolyPatch>(poly))
        {
            patchGIB[pI] = true;
            const indirectPolyPatch& inPoly =
                refCast<const indirectPolyPatch>(this->boundary()[pI].patch());
            const word& inPolyType =  inPoly.indirectPolyPatchType();
            if (inPolyType=="master")
            {
                masterGIB_ = pI;
            }
            else if (inPolyType=="slave")
            {
                slaveGIB_ = pI;
            }
        }
    }
    if (masterGIB_ < 0 || slaveGIB_ < 0)
    {
        FatalErrorInFunction
            << "GIB master and slave patches were not found in mesh."
            << nl << exit(FatalError);
    }
}


void Foam::dynamicGIBFvMesh::findGIBPatches(const word& zoneName)
{
    deleteDemandDrivenData(patchGIBPtr_);
    patchGIBPtr_ = new boolList(this->boundary().size(), false);
    boolList& patchGIB = *patchGIBPtr_;
    forAll(this->boundary(), pI)
    {
        const polyPatch& poly = this->boundary()[pI].patch();
        if (isA<indirectPolyPatch>(poly))
        {
            const indirectPolyPatch& inPoly =
                refCast<const indirectPolyPatch>(this->boundary()[pI].patch());

            const label& zoneId = inPoly.zoneId();
            if (this->faceZones()[zoneId].name() == zoneName)
            {
                patchGIB[pI] = true;
                const word& inPolyType =  inPoly.indirectPolyPatchType();
                if (inPolyType=="master")
                {
                    masterGIB_ = pI;
                }
                else if (inPolyType=="slave")
                {
                    slaveGIB_ = pI;
                }
            }
        }
    }
}

void Foam::dynamicGIBFvMesh::modifyRegionLabels(labelList& cellIndi) const
{
    const label fFluid = 0;
    const label fSolid = 1;
    labelList cellIndiBU = cellIndi;

    DynamicList<label> patchRegionLabel(this->boundary().size());

    forAll(region0Patch_, eI)
    {
        word substring = region0Patch_[eI];

        forAll(this->boundary(), patchI)
        {
            const fvPatch& cPatch(this->boundary()[patchI]);
            word patchName(cPatch.name());

            if (patchName.find(substring, 0) != string::npos)
            {
                const labelList& inCells = this->boundary()[patchI].faceCells();
                if (inCells.size()!=0)
                {
                    patchRegionLabel.append(cellIndiBU[inCells[0]]);
                }
            }
        }
    }
    patchRegionLabel.shrink();

    labelListList patchPro(Pstream::nProcs(), patchRegionLabel);

    Pstream::allGatherList(patchPro);

    forAll(cellIndi, cI)
    {
        bool found = false;
        forAll(patchPro, procI)
        {
            const labelList labels = patchPro[procI];
            forAll(labels, lI)
            {
                if (labels[lI] == cellIndiBU[cI])
                {
                    found = true;
                }
            }
        }
        if (!found)
        {
            cellIndi[cI] = fSolid;
        }
        else
        {
            cellIndi[cI] = fFluid;
        }
    }

    checkRegion(cellIndi);
}


bool Foam::dynamicGIBFvMesh::isIdenticalInterface() const
{
    bool comparison = true;

    if (fl() != fl0())
    {
        comparison = false;
    }

    reduce(comparison, andOp<bool>());

    return comparison;
}


void Foam::dynamicGIBFvMesh::smoothInternalBasePoints
(
    pointField& basePatchPoints,
    const labelList& pointAdd
) const
{
    const pointField& basePoints = *basePoints_;
    //-store the initial location of tge base interface in global addressing
    pointField gbasePatchPointsInit(basePoints);
    forAll(basePatchPoints, pI)
    {
        gbasePatchPointsInit[pointAdd[pI]] = basePatchPoints[pI];
    }
    syncPoints(gbasePatchPointsInit);

    //- mark points which are on faces connected on the boundary
    const labelList& interFaces = fl();
    const boolList& bPoints = boundaryPoints();

    //-smooth base patch with relaxation
    for (int i = 1; i<4; i++)
    {
        pointField gbasePatchPoints(basePoints);
        forAll(basePatchPoints, pI)
        {
            gbasePatchPoints[pointAdd[pI]] = basePatchPoints[pI];
        }

        //- construct patch
        indirectPrimitivePatch smoothPatch
        (
            IndirectList<face>(faces(), interFaces),
            gbasePatchPoints
        );

        PrimitivePatchInterpolation<indirectPrimitivePatch> pInterC
        (
            smoothPatch
        );

        tmp<pointField> smoothFaceAveraget =
            pInterC.faceToPointInterpolate(smoothPatch.faceCentres());
        pointField& smoothFaceAverage = smoothFaceAveraget.ref();

        const labelList& globalAddressing = smoothPatch.meshPoints();
        forAll(smoothFaceAverage, pI)
        {
            const label gpI = globalAddressing[pI];
            if (bPoints[gpI])
            {
                smoothFaceAverage[pI] = gbasePatchPoints[gpI];
            }
        }

        pointField gsmoothFaceAverage(basePoints);
        forAll(smoothFaceAverage, pI)
        {
            gsmoothFaceAverage[globalAddressing[pI]] = smoothFaceAverage[pI];
        }
        syncPoints(gsmoothFaceAverage);

        correctConstraintPatches(gsmoothFaceAverage);

        //- update the baseMesh points using relaxation
        forAll(basePatchPoints, pI)
        {
            const label gpI = pointAdd[pI];
            basePatchPoints[pI] =
                0.5*gbasePatchPointsInit[gpI]
              + 0.5*gsmoothFaceAverage[gpI];
        }
    }
}


void Foam::dynamicGIBFvMesh::updateSolidCellZone() const
{
    const labelList& cReg = cRegion();
    const cellZone& cZone = this->cellZones()[findOrCreateCellZone()];
    cellZone& ccZone =  const_cast<cellZone&>(cZone);

    DynamicList<label> zoneCL(cReg.size());

    forAll(cReg, cI)
    {
        if (cReg[cI]!=0)
        {
            zoneCL.append(cI);
        }
    }
    zoneCL.shrink();
    ccZone.clearAddressing();
    const labelList c(zoneCL);
    ccZone = c;
}


void Foam::dynamicGIBFvMesh::storeOldTimes()
{
    StoreOldTimeFields<scalar, fvPatchField, volMesh>(*this);
    StoreOldTimeFields<vector, fvPatchField, volMesh>(*this);
    StoreOldTimeFields<sphericalTensor, fvPatchField, volMesh>(*this);
    StoreOldTimeFields<symmTensor, fvPatchField, volMesh>(*this);
    StoreOldTimeFields<tensor, fvPatchField, volMesh>(*this);
    StoreOldTimeFields<scalar, fvsPatchField, surfaceMesh>(*this);
    StoreOldTimeFields<vector, fvsPatchField, surfaceMesh>(*this);
    StoreOldTimeFields<sphericalTensor, fvsPatchField, surfaceMesh>(*this);
    StoreOldTimeFields<symmTensor, fvsPatchField, surfaceMesh>(*this);
    StoreOldTimeFields<tensor, fvsPatchField, surfaceMesh>(*this);
}

void Foam::dynamicGIBFvMesh::correctBCs()
{
    // Only sync parallel BCs as there are issues correcting some BCs
    // at an intermediate stage
    CorrectProcessorBCs<scalar>(*this);
    CorrectProcessorBCs<vector>(*this);
    CorrectProcessorBCs<sphericalTensor>(*this);
    CorrectProcessorBCs<symmTensor>(*this);
    CorrectProcessorBCs<tensor>(*this);

    if (correctBCs_)
    {
/*
        CorrectBCs<scalar>(*this);
        CorrectBCs<vector>(*this);
        CorrectBCs<sphericalTensor>(*this);
        CorrectBCs<symmTensor>(*this);
        CorrectBCs<tensor>(*this);
*/
        if (this->foundObject<volScalarField>("k"))
        {
            volScalarField& k =
                const_cast<volScalarField&>
                (
                    lookupObject<volScalarField>("k")
                );
            forAll(k, cI)
            {
                if (k[cI] < SMALL)
                {
                    k[cI] = SMALL;
                }
            }
            volScalarField::Boundary& kbf = k.boundaryFieldRef();
            forAll(kbf, pI)
            {
                forAll(kbf[pI], pfI)
                {
                    if (kbf[pI][pfI] < SMALL)
                    {
                        kbf[pI][pfI] = SMALL;
                    }
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::explicitLimitingFields()
{
    if (limitNegFieldNames_.size())
    {
        forAll(limitNegFieldNames_, fieldI)
        {
            word substring = limitNegFieldNames_[fieldI];

            if (this->foundObject<volScalarField>(substring))
            {
                volScalarField& field = const_cast<volScalarField&>
                (
                    lookupObject<volScalarField>(substring)
                );
                label nT = field.nOldTimes();

                forAll(field, cI)
                {
                    if (field[cI]<0)
                    {
                        field[cI]=0;
                    }
                }
                if (nT!=0)
                {
                    forAll(field.oldTime(), cI)
                    {
                        if (field.oldTime()[cI]<0)
                        {
                            field.oldTime()[cI]=0;
                        }
                    }
                }
            }
            else
            {
                FatalErrorInFunction << abort(FatalError);
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::clearOutGIBData()
{
    deleteDemandDrivenData(multInterFacesPtr_);
    deleteDemandDrivenData(flPtr_);
    deleteDemandDrivenData(fl0Ptr_);
    deleteDemandDrivenData(fmPtr_);
    deleteDemandDrivenData(fm0Ptr_);
    deleteDemandDrivenData(interPointsPtr_);
    deleteDemandDrivenData(interPoints0Ptr_);
    deleteDemandDrivenData(phiBPtr_);
    deleteDemandDrivenData(phiTrPtr_);
    deleteDemandDrivenData(phiTrSPtr_);
    deleteDemandDrivenData(phiTrGPtr_);
    deleteDemandDrivenData(faceIndicatorPtr_);
    deleteDemandDrivenData(faceIndicator0Ptr_);
    deleteDemandDrivenData(cRegionPtr_);
    deleteDemandDrivenData(cRegion0Ptr_);
    deleteDemandDrivenData(recPoints0Ptr_);
    deleteDemandDrivenData(recAllPoints0Ptr_);
    deleteDemandDrivenData(hitPointPtr_);
    deleteDemandDrivenData(closeBoundaryPointsPtr_);
    deleteDemandDrivenData(markedBoundaryPointsPtr_);
    deleteDemandDrivenData(normal2DEdgesPtr_);
    deleteDemandDrivenData(boundaryPointsRevPtr_);
    deleteDemandDrivenData(fullSnapCellsPtr_);
    deleteDemandDrivenData(fullSnapCellPointsPtr_);
    deleteDemandDrivenData(popCellPointsPtr_);
    deleteDemandDrivenData(popSPointsPtr_);
    deleteDemandDrivenData(popGPointsPtr_);
    deleteDemandDrivenData(popUpCellsPtr_);
    deleteDemandDrivenData(meshQualityPointsPtr_);
    deleteDemandDrivenData(patchGIBPtr_);
    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(concavePointsPtr_);
    deleteDemandDrivenData(concaveEdgesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicGIBFvMesh::dynamicGIBFvMesh
(
    const IOobject& io,
    const word& typeN
)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dynamicMeshDict().subDict(typeN + "Coeffs")),
    basePoints_(nullptr),
    baseCC_(nullptr),
    baseSf_(nullptr),
    baseCf_(nullptr),
    triName_(dynamicMeshCoeffs_.lookup("triSurfaceName")),
    patchGIBPtr_(nullptr),
    masterGIB_(-1),
    slaveGIB_(-1),
    cRegionPtr_(nullptr),
    cRegion0Ptr_(nullptr),
    closeBoundaryPointsPtr_(nullptr),
    markedBoundaryPointsPtr_(nullptr),
    normal2DEdgesPtr_(nullptr),
    boundaryPointsRevPtr_(nullptr),
    ibMeshPtr_(nullptr),
    flPtr_(nullptr),
    fmPtr_(nullptr),
    interPointsPtr_(nullptr),
    interPoints0Ptr_(nullptr),
    boundaryPointsPtr_(nullptr),
    multInterFacesPtr_(nullptr),
    concavePointsPtr_(nullptr),
    concaveEdgesPtr_(nullptr),
    fullSnapCellsPtr_(nullptr),
    fullSnapCellPointsPtr_(nullptr),
    popCellPointsPtr_(nullptr),
    popSPointsPtr_(nullptr),
    popGPointsPtr_(nullptr),
    hitIndexPtr_(nullptr),
    hitPointPtr_(nullptr),
    popUpCellsPtr_(nullptr),
    faceIndicatorPtr_(nullptr),
    faceIndicator0Ptr_(nullptr),
    masterFCells_(boolList(this->cells().size(), false)),
    slaveFCells_(boolList(this->cells().size(), false)),
    region0Patch_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("region0Patch", wordList(1,"inlet"))
    ),
    limitNegFieldNames_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("limitNegFieldNames", wordList())
    ),
    popSubSteps_
    (
        dynamicMeshCoeffs_.lookupOrDefault<label>
        ("popSubSteps", 10)
    ),
    maxNonOrth_
    (
        cos
        (
            degToRad
            (
                dynamicMeshCoeffs_.lookupOrDefault<scalar>
                ("maxNonOrthogonality", 80)
            )
        )
    ),
    unsnapVar_
    (
        dynamicMeshCoeffs_.lookupOrDefault<scalar>
        ("unsnapCellsPercentage", 0.3)
    ),
    fl0Ptr_(nullptr),
    fm0Ptr_(nullptr),
    phiBPtr_(nullptr),
    phiTrPtr_(nullptr),
    phiTrSPtr_(nullptr),
    phiTrGPtr_(nullptr),
    debugMode_(dynamicMeshCoeffs_.lookupOrDefault<bool>("debug", false)),
    correctBCs_(dynamicMeshCoeffs_.lookupOrDefault<bool>("correctBCs", true)),
    pyrPrismFlip_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("flipPyramidsAndPrisms", false)
    ),
    boundPopValues_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("boundPopValues", true)
    ),
    checkInterface_(false),
    oldPoints_(points()),
    prevPoints_(points()),
    recPoints0Ptr_(nullptr),
    recAllPoints0Ptr_(nullptr),
    meshQualityPointsPtr_(nullptr)
{
    initialize();
}


Foam::dynamicGIBFvMesh::dynamicGIBFvMesh
(
    const IOobject& io,
    const dictionary dict
)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dict),
    basePoints_(nullptr),
    baseCC_(nullptr),
    baseSf_(nullptr),
    baseCf_(nullptr),
    triName_(dynamicMeshCoeffs_.lookup("triSurfaceName")),
    patchGIBPtr_(nullptr),
    masterGIB_(-1),
    slaveGIB_(-1),
    cRegionPtr_(nullptr),
    cRegion0Ptr_(nullptr),
    closeBoundaryPointsPtr_(nullptr),
    markedBoundaryPointsPtr_(nullptr),
    normal2DEdgesPtr_(nullptr),
    boundaryPointsRevPtr_(nullptr),
    ibMeshPtr_(nullptr),
    flPtr_(nullptr),
    fmPtr_(nullptr),
    interPointsPtr_(nullptr),
    interPoints0Ptr_(nullptr),
    boundaryPointsPtr_(nullptr),
    multInterFacesPtr_(nullptr),
    concavePointsPtr_(nullptr),
    concaveEdgesPtr_(nullptr),
    fullSnapCellsPtr_(nullptr),
    fullSnapCellPointsPtr_(nullptr),
    popCellPointsPtr_(nullptr),
    popSPointsPtr_(nullptr),
    popGPointsPtr_(nullptr),
    hitIndexPtr_(nullptr),
    hitPointPtr_(nullptr),
    popUpCellsPtr_(nullptr),
    faceIndicatorPtr_(nullptr),
    faceIndicator0Ptr_(nullptr),
    masterFCells_(boolList(this->cells().size(), false)),
    slaveFCells_(boolList(this->cells().size(), false)),
    region0Patch_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("region0Patch", wordList(1,"inlet"))
    ),
    limitNegFieldNames_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("limitNegFieldNames", wordList())
    ),
    popSubSteps_
    (
        dynamicMeshCoeffs_.lookupOrDefault<label>
        ("popSubSteps", 10)
    ),
    maxNonOrth_
    (
        cos
        (
            degToRad
            (
                dynamicMeshCoeffs_.lookupOrDefault<scalar>
                ("maxNonOrthogonality", 80)
            )
        )
    ),
    unsnapVar_
    (
        dynamicMeshCoeffs_.lookupOrDefault<scalar>
        ("unsnapCellsPercentage", 0.3)
    ),
    fl0Ptr_(nullptr),
    fm0Ptr_(nullptr),
    phiBPtr_(nullptr),
    phiTrPtr_(nullptr),
    phiTrSPtr_(nullptr),
    phiTrGPtr_(nullptr),
    debugMode_(dynamicMeshCoeffs_.lookupOrDefault<bool>("debug", false)),
    correctBCs_(dynamicMeshCoeffs_.lookupOrDefault<bool>("correctBCs", true)),
    pyrPrismFlip_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("flipPyramidsAndPrisms", false)
    ),
    boundPopValues_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("boundPopValues", true)
    ),
    checkInterface_(false),
    oldPoints_(points()),
    prevPoints_(points()),
    recPoints0Ptr_(nullptr),
    recAllPoints0Ptr_(nullptr),
    meshQualityPointsPtr_(nullptr)
{
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicGIBFvMesh::~dynamicGIBFvMesh()
{
    deleteDemandDrivenData(basePoints_);
    deleteDemandDrivenData(baseCC_);
    deleteDemandDrivenData(baseSf_);
    deleteDemandDrivenData(baseCf_);
    deleteDemandDrivenData(ibMeshPtr_);
    deleteDemandDrivenData(hitIndexPtr_);
    this->clearOutGIBData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicGIBFvMesh::updateInit(const word& zoneName)
{
    findGIBPatches(zoneName);
    clearOutGIBData();
    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = this->faceZones()[zoneId];
    faceZone& fZone =  const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(0);
    fm0Ptr_ = new boolList(0);

    tmp<pointField> snapPt = findSnappedPoints(false);
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList&  fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    if (tangentialBoundaryMotion())
    {
        doTangentialBoundaryMotion(snapP);
    }
    correctConstraintPatches(snapP);

    fvMesh::moveGIBPoints(snapP);
    updateGIB();

    word cRegionName = "cRegion";
    if (zoneName != word::null)
    {
        cRegionName += "_"+zoneName;
    }
    writeScalarField(cRegionName, cRegion(), true);
}


void Foam::dynamicGIBFvMesh::doUpdate(mapGIB& mapCl, bool fromBase)
{
    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = this->faceZones()[zoneId];
    faceZone& fZone =  const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    tmp<pointField> snapPt = findSnappedPoints(fromBase);
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList&  fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    if (tangentialBoundaryMotion())
    {
        doTangentialBoundaryMotion(snapP);
    }

    correctConstraintPatches(snapP);

    fvMesh::moveGIBPoints(snapP);

    checkInterface_ = isIdenticalInterface();

    if (!checkInterface_)
    {
        updateGIB();

        mapCl.mapBcs();

        popShrinkFields();

        correctV0();

        popGrowFields();
    }
    else
    {
        storeGIBFields();
        correctV0();
    }

    oldPoints_ = recAllPoints0();

    resetMeshFluxes();

    writeProblematicCells();

    faceCellsVisDebug();

    correctBCs();

    explicitLimitingFields();

    this->topoChanging(true);

}


void Foam::dynamicGIBFvMesh::writeGeometry(const fileName& surfaceName)
{
    const fvPatch& mPatch = this->boundary()[masterId()];

    MeshedSurface<face> surface
    (
        mPatch.patch(),
        true
    );

    if (Pstream::master())
    {

        fileName globalCasePath
        (
            surfaceName.isAbsolute()
          ? surfaceName
          : (
                this->time().processorCase()
              ? this->time().rootPath()/this->time().globalCaseName()/surfaceName
              : this->time().path()/surfaceName
            )
        );
        globalCasePath.clean();

        if (!exists(globalCasePath.path()))
        {
            Info<< "Path to surface does not exist." << endl;
            Info<< "Creating path: "
                 << globalCasePath.path()
                 << endl;

            mkDir(globalCasePath.path());
        }

        surface.write(globalCasePath);
    }
}


const Foam::labelList& Foam::dynamicGIBFvMesh::fl() const
{
    if (!flPtr_)
    {
        makeFl();
    }

    return *flPtr_;
}


const Foam::labelList& Foam::dynamicGIBFvMesh::fl0() const
{
    return *fl0Ptr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::fm() const
{
    if (!fmPtr_)
    {
        makeFm();
    }

    return *fmPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::fm0() const
{
    return *fm0Ptr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::interP() const
{
    if (!interPointsPtr_)
    {
        makeInterPoints();
    }

    return *interPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::interP0() const
{
    if (!interPoints0Ptr_)
    {
        makeInterPoints0();
    }
    return *interPoints0Ptr_;
}


const Foam::surfaceScalarField& Foam::dynamicGIBFvMesh::phiB() const
{
    if (!phiBPtr_)
    {
        makePhiB();
    }

    return *phiBPtr_;
}


const Foam::surfaceScalarField& Foam::dynamicGIBFvMesh::phiTr() const
{
    if (!phiTrPtr_)
    {
        makePhiTr();
    }

    return *phiTrPtr_;
}


const Foam::surfaceScalarField& Foam::dynamicGIBFvMesh::phiTrS() const
{
    if (!phiTrSPtr_)
    {
        makePhiTrS();
    }

    return *phiTrSPtr_;
}


const Foam::surfaceScalarField& Foam::dynamicGIBFvMesh::phiTrG() const
{
    if (!phiTrGPtr_)
    {
        makePhiTrG();
    }

    return *phiTrGPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::faceIndicator() const
{
    if (!faceIndicatorPtr_)
    {
        makeFaceIndicator();
    }

    return *faceIndicatorPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::faceIndicator0() const
{
    if (!faceIndicator0Ptr_)
    {
        makeFaceIndicator0();
    }

    return *faceIndicator0Ptr_;
}


const Foam::labelList& Foam::dynamicGIBFvMesh::cRegion() const
{
    if (!cRegionPtr_)
    {
        makecRegion();
    }

    return *cRegionPtr_;
}


const Foam::labelList& Foam::dynamicGIBFvMesh::cRegion0() const
{
    if (!cRegion0Ptr_)
    {
        makecRegion0();
    }

    return *cRegion0Ptr_;
}


const Foam::labelList& Foam::dynamicGIBFvMesh::fullSnapCells() const
{
    if (!fullSnapCellsPtr_)
    {
        makeFullSnapCells();
    }

    return *fullSnapCellsPtr_;
}


const Foam::labelList& Foam::dynamicGIBFvMesh::fullSnapCellPoints() const
{
    if (!fullSnapCellPointsPtr_)
    {
        makeFullSnapCellPoints();
    }

    return *fullSnapCellPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::popCellPoints() const
{
    if (!popCellPointsPtr_)
    {
        makePopCellPoints();
    }

    return *popCellPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::popSPoints() const
{
    if (!popSPointsPtr_)
    {
        makePopSPoints();
    }

    return *popSPointsPtr_;
}

const Foam::boolList& Foam::dynamicGIBFvMesh::popGPoints() const
{
    if (!popGPointsPtr_)
    {
        makePopGPoints();
    }

    return *popGPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::concavePoints() const
{
    if (!concavePointsPtr_)
    {
        makeConcavePoints();
    }

    return *concavePointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::concaveEdges() const
{
    if (!concaveEdgesPtr_)
    {
        makeConcavePoints();
    }

    return *concaveEdgesPtr_;
}


const Foam::labelListList& Foam::dynamicGIBFvMesh::closeBoundaryPoints() const
{
    if (!closeBoundaryPointsPtr_)
    {
        makeCloseBoundaryPoints();
    }

    return *closeBoundaryPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::markedBoundaryPoints() const
{
    if (!markedBoundaryPointsPtr_)
    {
        makeMarkedBoundaryPoints();
    }

    return *markedBoundaryPointsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::boundaryPointsRev() const
{
    if (!boundaryPointsRevPtr_)
    {
        makeBoundaryPointsRev();
    }

    return *boundaryPointsRevPtr_;
}


const Foam::pointField& Foam::dynamicGIBFvMesh::recPoints0() const
{
    if (!recPoints0Ptr_)
    {
        makeRecPoints0();
    }

    return *recPoints0Ptr_;
}


const Foam::pointField& Foam::dynamicGIBFvMesh::recAllPoints0() const
{
    if (!recAllPoints0Ptr_)
    {
        makeRecAllPoints0();
    }

    return *recAllPoints0Ptr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::popUpCells() const
{
    if (!popUpCellsPtr_)
    {
        makePopUpCells();
    }

    return *popUpCellsPtr_;
}


const Foam::boolList& Foam::dynamicGIBFvMesh::boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        makeBoundaryPoints();
    }

    return *boundaryPointsPtr_;
}


// ************************************************************************* //
