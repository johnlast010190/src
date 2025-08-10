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
    (c) 2016-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "meshes/Identifiers/surface/surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "indexedOctree/treeDataFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapGIB, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapGIB::storeFields()
{
    mapSerial();
}


void Foam::mapGIB::mapSerial()
{
    const label& mId = mesh_.masterId();
    const label& sId = mesh_.slaveId();
    StoreOldFieldsToPatch<scalar, fvPatchField, volMesh>(gmfvScalars_, mId);
    StoreOldFieldsToPatch<scalar, fvsPatchField, surfaceMesh>(gmfvsScalars_, mId);
    StoreOldFieldsToPatch<scalar, fvPatchField, volMesh>(gsfvScalars_, sId);
    StoreOldFieldsToPatch<scalar, fvsPatchField, surfaceMesh>(gsfvsScalars_, sId);

    StoreOldFieldsToPatch<vector, fvPatchField, volMesh>(gmfvVectors_, mId);
    StoreOldFieldsToPatch<vector, fvsPatchField, surfaceMesh>(gmfvsVectors_, mId);
    StoreOldFieldsToPatch<vector, fvPatchField, volMesh>(gsfvVectors_, sId);
    StoreOldFieldsToPatch<vector, fvsPatchField, surfaceMesh>(gsfvsVectors_, sId);

    StoreOldFieldsToPatch<tensor, fvPatchField, volMesh>(gmfvTensors_, mId);
    StoreOldFieldsToPatch<tensor, fvsPatchField, surfaceMesh>(gmfvsTensors_, mId);
    StoreOldFieldsToPatch<tensor, fvPatchField, volMesh>(gsfvTensors_, sId);
    StoreOldFieldsToPatch<tensor, fvsPatchField, surfaceMesh>(gsfvsTensors_, sId);

    StoreOldFieldsToPatch<sphericalTensor, fvPatchField, volMesh>(gmfvSpTensors_, mId);
    StoreOldFieldsToPatch<sphericalTensor, fvsPatchField, surfaceMesh>(gmfvsSpTensors_, mId);
    StoreOldFieldsToPatch<sphericalTensor, fvPatchField, volMesh>(gsfvSpTensors_, sId);
    StoreOldFieldsToPatch<sphericalTensor, fvsPatchField, surfaceMesh>(gsfvsSpTensors_, sId);

    StoreOldFieldsToPatch<symmTensor, fvPatchField, volMesh>(gmfvSymmTensors_, mId);
    StoreOldFieldsToPatch<symmTensor, fvsPatchField, surfaceMesh>(gmfvsSymmTensors_, mId);
    StoreOldFieldsToPatch<symmTensor, fvPatchField, volMesh>(gsfvSymmTensors_, sId);
    StoreOldFieldsToPatch<symmTensor, fvsPatchField, surfaceMesh>(gsfvsSymmTensors_, sId);
}

void Foam::mapGIB::combinePolyPatch()
{
    sourcePointsm_ = sourcePatch_.localPoints();

    sourceFacesm_ = sourcePatch_.localFaces();

    const fvPatch& gibPatch(mesh_.boundary()[mesh_.slaveId()]);

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    faceList faces(gibPolyPatch);
    {
        const boolList& flipmap = gibPolyPatch.fm();
        forAll(faces, fI)
        {
            if (!flipmap[fI])
            {
                faces[fI].flip();
            }
        }
    }

    sourceFacess_ = faces;
    sourcePointss_ = gibPolyPatch.points();

    sourcePatchm_ = new primitivePatch
    (
        SubList<face>(sourceFacesm_, sourceFacesm_.size()),
        sourcePointsm_
    );

    sourcePatchs_ = new primitivePatch
    (
        SubList<face>(sourceFacess_, sourceFacess_.size()),
        sourcePointss_
    );
    if (false)
    {
        simpleVTKWriter writeMaster
        (
            sourcePatchm_->localFaces(),
            sourcePatchm_->localPoints()
        );
        writeMaster.addFaceData("normal", sourcePatchm_->faceNormals());

        simpleVTKWriter writeSlave
        (
            sourcePatchs_->localFaces(),
            sourcePatchs_->localPoints()
        );
        writeSlave.addFaceData("normal", sourcePatchs_->faceNormals());

        simpleVTKWriter writeSlave2
        (
            gibPatch.patch().localFaces(),
            gibPatch.patch().localPoints()
        );

        writeMaster.write("master.vtk");
        writeSlave.write("slave.vtk");
        writeSlave2.write("slave2.vtk");
    }


    if (Pstream::parRun())
    {
        List<Field<vector >> ppPoints(Pstream::nProcs());
        List<faceList> ppFaces(Pstream::nProcs());

        ppPoints[Pstream::myProcNo()] = sourcePatch_.localPoints();
        ppFaces[Pstream::myProcNo()] = sourcePatch_.localFaces();

        Pstream::allGatherList(ppPoints);
        Pstream::allGatherList(ppFaces);

        int iPro=1;
        int nPro=Pstream::nProcs();
        label nPoints = 0;
        for (int i = iPro; i < nPro; i++)
        {
            nPoints +=
            ppPoints[i-1].size();
            forAll(ppFaces[i], fI)
            {
                forAll(ppFaces[i][fI], pI)
                {
                    ppFaces[i][fI][pI] += nPoints;
                }
            }
        }

        gppPoints_ = ListListOps::combine<Field<vector>>
            (
                ppPoints,
                accessOp<Field<vector>>()
            );

        gppFaces_ = ListListOps::combine<faceList>
            (
                ppFaces, accessOp<faceList>()
            );
    }
    else
    {
        gppPoints_ = sourcePatch_.localPoints();
        gppFaces_ = sourcePatch_.localFaces();
    }
}


void Foam::mapGIB::mapMaster()
{
    const label& mId = mesh_.masterId();
    const indirectPolyPatch& masterGIB
    (
        refCast<const indirectPolyPatch>
        (
            mesh_.boundary()[mId].patch()
        )
    );
    const boolList& flipmap = masterGIB.fm();

    targetPointsm_ = masterGIB.localPoints();
    targetFacesm_ = masterGIB.localFaces();

    forAll(targetFacesm_, fI)
    {
        if (flipmap[fI])
        {
            targetFacesm_[fI].flip();
        }
    }

    targetPatchm_ = new primitivePatch
    (
        SubList<face>(targetFacesm_, targetFacesm_.size()),
        targetPointsm_
    );

    const AMIPatchToPatchInterpolation& gibmasterIn = masterAMIInter();

    if (false)
    {
        simpleVTKWriter writeMasterSource
        (
            sourcePatchm_->localFaces(),
            sourcePatchm_->localPoints()
        );
        writeMasterSource.addFaceData("normal", sourcePatchm_->faceNormals());
        simpleVTKWriter writeMasterTarget
        (
            targetPatchm_->localFaces(),
            targetPatchm_->localPoints()
        );
        writeMasterTarget.addFaceData("normal", targetPatchm_->faceNormals());

        writeMasterSource.write("source_master.vtk");
        writeMasterTarget.write("target_master.vtk");
    }

    MapGIBField<scalar, fvPatchField, volMesh>
    (
        gmfvScalars_, mId, gibmasterIn
    );
    MapGIBField<vector, fvPatchField, volMesh>
    (
        gmfvVectors_, mId, gibmasterIn
    );
    MapGIBField<sphericalTensor, fvPatchField, volMesh>
    (
        gmfvSpTensors_, mId, gibmasterIn
    );
    MapGIBField<symmTensor, fvPatchField, volMesh>
    (
        gmfvSymmTensors_, mId, gibmasterIn
    );
    MapGIBField<tensor, fvPatchField, volMesh>
    (
        gmfvTensors_, mId, gibmasterIn
    );

    MapGIBField<scalar, fvsPatchField, surfaceMesh>
    (
        gmfvsScalars_, mId, gibmasterIn
    );
    MapGIBField<vector, fvsPatchField, surfaceMesh>
    (
        gmfvsVectors_, mId, gibmasterIn
    );
    MapGIBField<sphericalTensor, fvsPatchField, surfaceMesh>
    (
        gmfvsSpTensors_, mId, gibmasterIn
    );
    MapGIBField<symmTensor, fvsPatchField, surfaceMesh>
    (
        gmfvsSymmTensors_, mId, gibmasterIn
    );
    MapGIBField<tensor, fvsPatchField, surfaceMesh>
    (
        gmfvsTensors_, mId, gibmasterIn
    );
}


void Foam::mapGIB::mapSlave()
{
    const label& sId = mesh_.slaveId();
    const indirectPolyPatch& slaveGIB
    (
        refCast<const indirectPolyPatch>
        (
            mesh_.boundary()[sId].patch()
        )
    );
    const boolList& flipmap = slaveGIB.fm();

    targetPointss_ = slaveGIB.localPoints();
    targetFacess_ = slaveGIB.localFaces();

    forAll(targetFacess_, fI)
    {
        if (!flipmap[fI])
        {
            targetFacess_[fI].flip();
        }
    }

    targetPatchs_ = new primitivePatch
    (
        SubList<face>(targetFacess_, targetFacess_.size()),
        targetPointss_
    );

    const AMIPatchToPatchInterpolation& gibslaveIn = slaveAMIInter();

    if (false)
    {
        simpleVTKWriter writeSlaveSource
        (
            sourcePatchs_->localFaces(),
            sourcePatchs_->localPoints()
        );
        writeSlaveSource.addFaceData("normal", sourcePatchs_->faceNormals());
        simpleVTKWriter writeSlaveTarget
        (
            targetPatchs_->localFaces(),
            targetPatchs_->localPoints()
        );
        writeSlaveTarget.addFaceData("normal", targetPatchs_->faceNormals());

        writeSlaveSource.write("source_slave.vtk");
        writeSlaveTarget.write("target_slave.vtk");
    }

    MapGIBField<scalar, fvPatchField, volMesh>
    (
        gsfvScalars_, sId, gibslaveIn
    );
    MapGIBField<vector, fvPatchField, volMesh>
    (
        gsfvVectors_, sId, gibslaveIn
    );
    MapGIBField<sphericalTensor, fvPatchField, volMesh>
    (
        gsfvSpTensors_, sId, gibslaveIn
    );
    MapGIBField<symmTensor, fvPatchField, volMesh>
    (
        gsfvSymmTensors_, sId, gibslaveIn
    );
    MapGIBField<tensor, fvPatchField, volMesh>
    (
        gsfvTensors_, sId, gibslaveIn
    );

    MapGIBField<scalar, fvsPatchField, surfaceMesh>
    (
        gsfvsScalars_, sId, gibslaveIn
    );
    MapGIBField<vector, fvsPatchField, surfaceMesh>
    (
        gsfvsVectors_, sId, gibslaveIn
    );
    MapGIBField<sphericalTensor, fvsPatchField, surfaceMesh>
    (
        gsfvsSpTensors_, sId, gibslaveIn
    );
    MapGIBField<symmTensor, fvsPatchField, surfaceMesh>
    (
        gsfvsSymmTensors_, sId, gibslaveIn
    );
    MapGIBField<tensor, fvsPatchField, surfaceMesh>
    (
        gsfvsTensors_, sId, gibslaveIn
    );
}


void Foam::mapGIB::mapNonOverlappingTargetFaces()
{
    {
        const scalarField& tgtSum = masterAMIInter().tgtWeightsSum();
        DynamicList<label> dNonOverFaces(tgtSum.size());

        forAll(tgtSum, fI)
        {
            if (tgtSum[fI] == 0)
            {
                dNonOverFaces.append(fI);
            }
        }

        dNonOverFaces.shrink();

        // All processors need to take same path due to a gAverage in
        // MapNonOverlapFaces
        label nNonOverFaces = dNonOverFaces.size();
        reduce(nNonOverFaces, sumOp<label>());

        if (nNonOverFaces)
        {
            const primitivePatch& targetPm = *targetPatchm_;
            const primitivePatch& sourcePm = *sourcePatchm_;
            const label& mId = mesh_.masterId();
            const indirectPolyPatch& masterGIB
            (
                refCast<const indirectPolyPatch>
                (
                    mesh_.boundary()[mId].patch()
                )
            );
            const labelList patchFaces = masterGIB.fAddr();

            labelList indexes  = calcNearestPatchFaceMapping
                (
                    dNonOverFaces,
                    targetPm,
                    sourcePm,
                    patchFaces
                );

            mapNonOverlapFaces(dNonOverFaces, indexes, mId);
        }
    }

    {
        const scalarField& tgtSum = slaveAMIInter().tgtWeightsSum();
        DynamicList<label> dNonOverFaces(tgtSum.size());

        forAll(tgtSum, fI)
        {
            if (tgtSum[fI] == 0)
            {
                dNonOverFaces.append(fI);
            }
        }

        dNonOverFaces.shrink();

        // All processors need to take same path due to a gAverage in
        // MapNonOverlapFaces
        label nNonOverFaces = dNonOverFaces.size();
        reduce(nNonOverFaces, sumOp<label>());

        if (nNonOverFaces)
        {
            const primitivePatch& targetPs = *targetPatchs_;
            const primitivePatch& sourcePs = *sourcePatchs_;
            const label& sId = mesh_.slaveId();
            const indirectPolyPatch& slaveGIB
            (
                refCast<const indirectPolyPatch>
                (
                    mesh_.boundary()[sId].patch()
                )
            );
            const labelList patchFaces = slaveGIB.fAddr();

            labelList indexes  = calcNearestPatchFaceMapping
                (
                    dNonOverFaces,
                    targetPs,
                    sourcePs,
                    patchFaces
                );

            mapNonOverlapFaces(dNonOverFaces, indexes, sId);
        }
    }
}


void Foam::mapGIB::makeMasterAMIPatchToPatchInterpolation() const
{
    if (mAMIInterPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    mAMIInterPtr_ = new AMIPatchToPatchInterpolation
        (
            *sourcePatchm_,
            *targetPatchm_,
            faceAreaIntersect::tmMesh,
            false,
            //"partialFaceAreaWeightAMI",
            "faceAreaWeightAMI",
            -1,
            true,
            debug::floatOptimisationSwitch("AMIRelTol", 0.01),
            0.0,
            AMICosMatchAngle_,
            maxAMIWeightScale_
        );

}


Foam::labelList Foam::mapGIB::calcNearestPatchFaceMapping
(
    const labelList& dNonOverFaces,
    const primitivePatch& targetPm,
    const primitivePatch& sourcePm,
    const labelList& patchFaces
)const
{
    typedef Tuple2<pointIndexHit, Tuple2<scalar, label>> nearInfo;

    Random rndGen(123456);

    treeBoundBox patchBb
    (
        treeBoundBox(sourcePm.points(), sourcePm.meshPoints()).extend
        (
            rndGen,
            1e-4
        )
    );
    patchBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    patchBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    indexedOctree<treeDataPrimitivePatch<primitivePatch>> boundaryTree
    (
        treeDataPrimitivePatch<primitivePatch>
        (
            false,      // do not cache bb
            sourcePm,
            indexedOctree<treeDataPrimitivePatch<primitivePatch>>::perturbTol()
        ),
        patchBb,        // overall search domain
        8,              // maxLevel
        10,             // leafsize
        3.0             // duplicity
    );

    List<nearInfo> nearest(dNonOverFaces.size());

    labelList indexes= labelList(dNonOverFaces.size(), -1);

    forAll(dNonOverFaces, fI)
    {
        label fII = dNonOverFaces[fI];
        const point& sample = targetPm.faceCentres()[fII];

        pointIndexHit& nearInfo = nearest[fI].first();
        nearInfo = boundaryTree.findNearest
        (
            sample,
            magSqr(patchBb.span())
        );

        if (nearInfo.hit())
        {
            point fc(sourcePm[nearInfo.index()].centre(sourcePm.points()));
            nearInfo.setPoint(fc);
            nearest[fI].second().first() = magSqr(fc-sample);
            nearest[fI].second().second() = Pstream::myProcNo();
            indexes[fI] = nearInfo.index();
        }
        else
        {
            if (sourcePm.size()==0)
            {
                WarningInFunction << "Source patch size for mapping zero"
                                  << endl;

            }
            WarningInFunction << "Unable to map face " << fII << " at "
                              << targetPatchm_->faceCentres()[fII]
                              << endl;
        }
    }
    if (debug)
    {
        /*
        Info<< "Checking mapping weights" <<endl;
        Info<< "faceI" << tab
             << "weightSum" << tab
             << "faceCentre" << tab
             << endl;
        */


        forAll(dNonOverFaces, fI)
        {
            label fII = dNonOverFaces[fI];
            Pout<< "faceI" << tab
                 << "weightSum" << tab
                 << "faceCentre" << tab
                 << endl;
            Pout<< "target" << tab
                 << fII << tab
                 //<< tgtSum[fII] << tab
                 << targetPatchm_->faceCentres()[fII] << tab
                 << endl;
        }
    }
    return indexes;
}


void Foam::mapGIB::makeSlaveAMIPatchToPatchInterpolation() const
{
    if (sAMIInterPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    sAMIInterPtr_ = new AMIPatchToPatchInterpolation
        (
            *sourcePatchs_,
            *targetPatchs_,
            faceAreaIntersect::tmMesh,
            false,
            //"partialFaceAreaWeightAMI",
            "faceAreaWeightAMI",
            -1,
            true,
            debug::floatOptimisationSwitch("AMIRelTol", 0.01),
            0.0,
            AMICosMatchAngle_,
            maxAMIWeightScale_
        );
}


void Foam::mapGIB::clearOutGIBData()
{
    deleteDemandDrivenData(mAMIInterPtr_);
    deleteDemandDrivenData(sAMIInterPtr_);
    deleteDemandDrivenData(sourcePatchm_);
    deleteDemandDrivenData(sourcePatchs_);
    deleteDemandDrivenData(targetPatchm_);
    deleteDemandDrivenData(targetPatchs_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapGIB::mapGIB
(
    const dynamicGIBFvMesh& mesh,
    const primitivePatch& sourcePatch,
    bool mapSlave,
    const scalar& AMICosMatchAngle,
    const scalar& maxAMIWeightScale
)
:
    mesh_(mesh),
    mapSlave_(mapSlave),
    AMICosMatchAngle_(AMICosMatchAngle),
    maxAMIWeightScale_(maxAMIWeightScale),
    gmfvScalars_(0),
    gmfvVectors_(0),
    gmfvSpTensors_(0),
    gmfvSymmTensors_(0),
    gmfvTensors_(0),
    gmfvsScalars_(0),
    gmfvsVectors_(0),
    gmfvsSpTensors_(0),
    gmfvsSymmTensors_(0),
    gmfvsTensors_(0),
    gsfvScalars_(0),
    gsfvVectors_(0),
    gsfvSpTensors_(0),
    gsfvSymmTensors_(0),
    gsfvTensors_(0),
    gsfvsScalars_(0),
    gsfvsVectors_(0),
    gsfvsSpTensors_(0),
    gsfvsSymmTensors_(0),
    gsfvsTensors_(0),
    mAMIInterPtr_(nullptr),
    sAMIInterPtr_(nullptr),
    sourcePatch_
    (
        sourcePatch
    ),
    sourcePatchm_(nullptr),
    sourcePatchs_(nullptr),
    sourcePointsm_(0),
    sourcePointss_(0),
    sourceFacesm_(0),
    sourceFacess_(0),
    targetPatchm_(nullptr),
    targetPatchs_(nullptr),
    targetPointsm_(0),
    targetPointss_(0),
    targetFacesm_(0),
    targetFacess_(0),
    gppPoints_(0),
    gppFaces_(0)
{
    combinePolyPatch();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mapGIB::mapBcs()
{
    storeFields();

    mapMaster();

    if (mapSlave_)
    {
        mapSlave();
    }
    mapNonOverlappingTargetFaces();
}

void Foam::mapGIB::map1To1()
{
    storeFields();

    Map1to1GIBField<scalar, fvPatchField, volMesh>
    (
        gmfvScalars_, gsfvScalars_
    );
    Map1to1GIBField<vector, fvPatchField, volMesh>
    (
        gmfvVectors_, gsfvVectors_
    );
    Map1to1GIBField<sphericalTensor, fvPatchField, volMesh>
    (
        gmfvSpTensors_, gsfvSpTensors_
    );
    Map1to1GIBField<symmTensor, fvPatchField, volMesh>
    (
        gmfvSymmTensors_, gsfvSymmTensors_
    );
    Map1to1GIBField<tensor, fvPatchField, volMesh>
    (
        gmfvTensors_, gsfvTensors_
    );

    Map1to1GIBField<scalar, fvsPatchField, surfaceMesh>
    (
        gmfvsScalars_, gsfvsScalars_
    );
    Map1to1GIBField<vector, fvsPatchField, surfaceMesh>
    (
        gmfvsVectors_, gsfvsVectors_
    );
    Map1to1GIBField<sphericalTensor, fvsPatchField, surfaceMesh>
    (
        gmfvsSpTensors_, gsfvsSpTensors_
    );
    Map1to1GIBField<symmTensor, fvsPatchField, surfaceMesh>
    (
        gmfvsSymmTensors_, gsfvsSymmTensors_
    );
    Map1to1GIBField<tensor, fvsPatchField, surfaceMesh>
    (
        gmfvsTensors_, gsfvsTensors_
    );
}

void Foam::mapGIB::mapNonOverlapFaces
(
    const labelList& dNonOverFaces,
    const labelList& nonOverFaces,
    const label patchID
)
{
    const label& mId = mesh_.masterId();
    const label& sId = mesh_.slaveId();
    if (mId == patchID)
    {
        MapNonOverlapFaces<scalar, fvPatchField, volMesh>
        (
            gmfvScalars_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<vector, fvPatchField, volMesh>
        (
            gmfvVectors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<sphericalTensor, fvPatchField, volMesh>
        (
            gmfvSpTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<symmTensor, fvPatchField, volMesh>
        (
            gmfvSymmTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<tensor, fvPatchField, volMesh>
        (
            gmfvTensors_, patchID, dNonOverFaces, nonOverFaces
        );

        MapNonOverlapFaces<scalar, fvsPatchField, surfaceMesh>
        (
            gmfvsScalars_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<vector, fvsPatchField, surfaceMesh>
        (
            gmfvsVectors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<sphericalTensor, fvsPatchField, surfaceMesh>
        (
            gmfvsSpTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<symmTensor, fvsPatchField, surfaceMesh>
        (
            gmfvsSymmTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<tensor, fvsPatchField, surfaceMesh>
        (
            gmfvsTensors_, patchID, dNonOverFaces, nonOverFaces
        );
    }
    else if (sId == patchID)
    {
        MapNonOverlapFaces<scalar, fvPatchField, volMesh>
        (
            gsfvScalars_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<vector, fvPatchField, volMesh>
        (
            gsfvVectors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<sphericalTensor, fvPatchField, volMesh>
        (
            gsfvSpTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<symmTensor, fvPatchField, volMesh>
        (
            gsfvSymmTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<tensor, fvPatchField, volMesh>
        (
            gsfvTensors_, patchID, dNonOverFaces, nonOverFaces
        );

        MapNonOverlapFaces<scalar, fvsPatchField, surfaceMesh>
        (
            gsfvsScalars_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<vector, fvsPatchField, surfaceMesh>
        (
            gsfvsVectors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<sphericalTensor, fvsPatchField, surfaceMesh>
        (
            gsfvsSpTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<symmTensor, fvsPatchField, surfaceMesh>
        (
            gsfvsSymmTensors_, patchID, dNonOverFaces, nonOverFaces
        );
        MapNonOverlapFaces<tensor, fvsPatchField, surfaceMesh>
        (
            gsfvsTensors_, patchID, dNonOverFaces, nonOverFaces
        );
    }
}

void Foam::mapGIB::mapInterfaceVelocity(scalarField* pG)
{
/*
    const patchToIndirectPatchInterpolation& masterIn = masterInter();
    tmp<scalarField > pGt =
        masterIn.pointInterpolate(*pG);

    const label& mId = mesh_.masterId();
    const labelList& mPoints = mesh_.boundary()[mId].patch().meshPoints();
    scalarField Gglobal = scalarField(mesh_.points().size() ,0);
    forAll(mPoints, pI)
    {
        Gglobal[mPoints[pI]] = pGt()[pI];
    }

    forAll(mPoints, pI)
    {
        pGt()[pI] = Gglobal[mPoints[pI]];
    }

    delete pG;
    pG = new scalarField(pGt());
    */
}


Foam::triSurface Foam::mapGIB::triS()
{
    pointField gibPoints = gppPoints_;
    faceList gibfaces = gppFaces_;

    labelList allZones(gppFaces_.size(), mesh_.slaveId());
    surfZoneIdentifierList surfZones(1);

    UnsortedMeshedSurface<face> unsortedFace
    (
        xferMove(gibPoints),
        xferMove(gibfaces),
        xferMove(allZones),
        xferMove(surfZones)
    );
    MeshedSurface<face> sortedFace(unsortedFace);

    labelList triMap;
    sortedFace.triangulate(triMap);


    triFaceList tFaces (sortedFace.surfFaces().size());
    forAll(tFaces, tFI)
    {
        point pA = sortedFace.points()[sortedFace.surfFaces()[tFI][0]];
        point pB = sortedFace.points()[sortedFace.surfFaces()[tFI][1]];
        point pC = sortedFace.points()[sortedFace.surfFaces()[tFI][2]];
        vector vAB = pB-pA;
        vector vBC = pC-pB;
        vector vABC = vAB^vBC;
        //scalar dotSf = gibPatch.Sf()[triMap[tFI]]&vABC;
        scalar dotSf = vector(0,1,0)&vABC;
        if (dotSf>0)
        {

            tFaces[tFI] = triFace
            (
                 sortedFace.surfFaces()[tFI][0],
                 sortedFace.surfFaces()[tFI][2],
                 sortedFace.surfFaces()[tFI][1]
            );

        }
        else
        {

            tFaces[tFI] = triFace
            (
                 sortedFace.surfFaces()[tFI][0],
                 sortedFace.surfFaces()[tFI][1],
                 sortedFace.surfFaces()[tFI][2]
            );

        }
    }
    triSurface newTriS = triSurface(tFaces, sortedFace.points());
    orientedSurface orSurf(newTriS, false);
    orSurf.orientConsistent(newTriS);
//    newTriS.cleanup(true);
    return newTriS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
