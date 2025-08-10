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
    (c) 2020 ESI

\*---------------------------------------------------------------------------*/

#include "foamVDB.H"
#include <openvdb/tools/ValueTransformer.h> //for foreach
#include <tbb/parallel_scan.h>


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// TBB body object for threaded copy of OpenFOAM List to VDB vector
template<class FoamType, class VDBType>
struct FoamToVDBOp
{
    const Foam::List<FoamType>& listIn_;
    std::vector<VDBType>& listOut_;

    FoamToVDBOp
    (
        const Foam::List<FoamType>& listIn,
        std::vector<VDBType>& listOut
    )
    :
        listIn_(listIn),
        listOut_(listOut)
    {
    }

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t n = range.begin(); n < range.end(); ++n)
        {
            VDBType outValue(VDBType::zero());

            for (size_t i = 0; i < VDBType::size; i++)
            {
                outValue[i] = listIn_[n][i];
            }

            listOut_[n] = outValue;
        }
    }
}; //FoamToVDBOp


// TBB body object for threaded send grids to neighour nodes
template<typename GridType>
struct SendGridOp
{
    const typename GridType::Ptr grid_;
    Foam::PstreamBuffers& pBufSize_;
    Foam::PstreamBuffers& pBufGrid_;
    const Foam::boundBox& myProcBB_;
    const Foam::List<Foam::boundBox>& procsBB_;
    bool isCellLevelGrid_;

    SendGridOp
    (
        const typename GridType::Ptr grid,
        Foam::PstreamBuffers& pBufSize,
        Foam::PstreamBuffers& pBufGrid,
        const Foam::boundBox& myProcBB,
        const Foam::List<Foam::boundBox>& procsBB,
        bool isCellLevelGrid
    )
    :
        grid_(grid),
        pBufSize_(pBufSize),
        pBufGrid_(pBufGrid),
        myProcBB_(myProcBB),
        procsBB_(procsBB),
        isCellLevelGrid_(isCellLevelGrid)
    {}

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t proci = range.begin(); proci < range.end(); ++proci)
        {
            if (proci == size_t(Foam::Pstream::myProcNo())) continue;

            const Foam::boundBox& otherProcBB = procsBB_[proci];

            if
            (
                myProcBB_.overlaps(otherProcBB)
             && (
                   isCellLevelGrid_
                 ? proci > size_t(Foam::Pstream::myProcNo())
                 : true
                )
            )
            {
                Foam::foamVDB::sendGrid<GridType>
                (
                    grid_,
                    pBufSize_,
                    pBufGrid_,
                    proci
                );
            } //if overlaps
        }
    }
}; // SendGridOp


// TBB body object for threaded receive grids from neighour nodes
template<typename GridType>
struct ReceiveGridsOp
{
    std::vector<typename GridType::Ptr>& gridsFromProcs_;
    Foam::PstreamBuffers& pBufSize_;
    Foam::PstreamBuffers& pBufGrid_;
    const Foam::boundBox& myProcBB_;
    const Foam::List<Foam::boundBox>& procsBB_;
    bool isCellLevelGrid_;

    ReceiveGridsOp
    (
        std::vector<typename GridType::Ptr>& gridsFromProcs,
        Foam::PstreamBuffers& pBufSize,
        Foam::PstreamBuffers& pBufGrid,
        const Foam::boundBox& myProcBB,
        const Foam::List<Foam::boundBox>& procsBB,
        bool isCellLevelGrid
    )
    :
        gridsFromProcs_(gridsFromProcs),
        pBufSize_(pBufSize),
        pBufGrid_(pBufGrid),
        myProcBB_(myProcBB),
        procsBB_(procsBB),
        isCellLevelGrid_(isCellLevelGrid)
    {}

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t proci = range.begin(); proci < range.end(); ++proci)
        {
            size_t myProcNo = Foam::Pstream::myProcNo();

            if (proci == myProcNo) continue;

            const Foam::boundBox& otherProcBB = procsBB_[proci];

            if
            (
                myProcBB_.overlaps(otherProcBB)
             && (
                   isCellLevelGrid_
                 ? proci < myProcNo
                 : true
                )
            )
            {
                typename GridType::Ptr gridFromProcI =
                    Foam::foamVDB::receiveGrid<GridType>
                    (
                        pBufSize_,
                        pBufGrid_,
                        proci
                    );

                gridsFromProcs_[proci] = gridFromProcI;
            } //if overlaps
        }
    }
}; // ReceiveGridsOp


// TBB body object for threaded receive grids from neighour nodes
template<typename GridType>
struct GatherGridsOp
{
    std::vector<typename GridType::Ptr>& gridsFromProcs_;
    Foam::PstreamBuffers& pBufSize_;
    Foam::PstreamBuffers& pBufGrid_;

    GatherGridsOp
    (
        std::vector<typename GridType::Ptr>& gridsFromProcs,
        Foam::PstreamBuffers& pBufSize,
        Foam::PstreamBuffers& pBufGrid
    )
    :
        gridsFromProcs_(gridsFromProcs),
        pBufSize_(pBufSize),
        pBufGrid_(pBufGrid)
    {}

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t proci = range.begin(); proci < range.end(); ++proci)
        {
            if (proci == Foam::Pstream::myProcNo()) continue;

            typename GridType::Ptr gridFromProcI =
                Foam::foamVDB::receiveGrid<GridType>
                (
                    pBufSize_,
                    pBufGrid_,
                    proci
                );

            gridsFromProcs_[proci] = gridFromProcI;
        }
    }
}; // GatherGridsOp


template<class FoamType, class VDBType>
void Foam::foamVDB::foamToVDB
(
    const List<FoamType>& inputList,
    std::vector<VDBType>& outputList
)
{
    if (VDBType::size != inputList[0].size())
    {
        FatalErrorInFunction
            << "OpenVDB data size (" << VDBType::size << ")"
            << "different from field size (" << inputList[0].size() << ")"
            << exit(FatalError);
    }

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, inputList.size()),
        FoamToVDBOp<FoamType, VDBType>
        (
            inputList,
            outputList
        )
    );
} // foamToVDB


template<class T>
Foam::label Foam::foamVDB::activeCount
(
    std::vector<T>& list
)
{
    const auto grainSize =
        std::max<size_t>
        (
            list.size() / tbb::this_task_arena::max_concurrency(),
            1024
        );
    const tbb::blocked_range<size_t> range(0, list.size(), grainSize);

    label nInactive =
        tbb::parallel_scan
        (
            range,
            /*identity*/0,
            /*scan body*/
            [&](const tbb::blocked_range<size_t>& r, const label& init, bool isFinalScan)
            {
                label res = init;
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    if (!list[i].active_)
                    {
                        res++;
                    }

                    if (isFinalScan && list[i].active_)
                    {
                        list[i].id_ -= res;
                    }
                }
                return res;
            },
            /*combine body*/
            [](const label& x, const label& y)
            {
                return x+y;
            },
            tbb::simple_partitioner()
        );

    //update id of inactive points
    if (typeid(T) == typeid(vdbPoint))
    {
        tbb::parallel_for
        (
            range,
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    if (!list[i].active_)
                    {
                        label oldID = list[i].id_;
                        list[i].id_ = list[oldID].id_;
                    }
                }
            },
            tbb::simple_partitioner()
        );
    }

    return list.size() - nInactive;
} //activeCount

template<typename GridType>
bool Foam::foamVDB::isValid
(
    typename GridType::ConstPtr grid
)
{
    if (!grid) return false;

    if (grid->activeVoxelCount() == 0) return false;

    return true;
}


template<typename GridType>
void Foam::foamVDB::trimGrids
(
    const std::vector<typename GridType::Ptr>& grids
)
{
    Timer timer("trimGrids");

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                if (!grids[cellLevel])  continue;

                if (grids[cellLevel+1])
                {
                    //TODO parallel_for on leafNode range
                    openvdb::tools::foreach
                    (
                        grids[cellLevel]->cbeginValueOn(),
                        [&](const typename GridType::ValueOnCIter& iter)
                        {
                            typename GridType::Accessor fineGridAcc
                            (
                                grids[cellLevel+1]->getAccessor()
                            );

                            const openvdb::Coord& ijk = iter.getCoord();

                            for (size_t i = 0; i < 8; i++)
                            {
                                const openvdb::Coord ijkFine = (ijk << 1) + COARSE_TO_FINE[i];

                                fineGridAcc.setValueOff(ijkFine);
                            }
                        },
                        /*threaded*/true,
                        /*shareOp*/false
                    );
                }
            } //for cellLevel
        } //TBB body trim
    );
} //trimGrids

template<typename GridType>
void Foam::foamVDB::sendGrids
(
    const std::vector<typename GridType::Ptr>& inputGrids,
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    std::ostringstream ostr(std::ios_base::binary);

    openvdb::GridPtrVecPtr grids(new openvdb::GridPtrVec);

    for (size_t i = 0; i < inputGrids.size(); i++)
    {
        grids->push_back(inputGrids[i]);
    }

    openvdb::io::Stream(ostr).write(*grids);

    std::streamsize strSize = ostr.str().size();

    std::string gridName = inputGrids[0]->getName();

    //Pout<< "[thread ID " << Foam::name(tbb::task_arena::current_thread_index()) << "] "
    //    << "Sending " << gridName
    //    << " to node " << node
    //    << " - size " << strSize/1000 << " kB"
    //    << endl;

    UOPstream toProcSize(node, pBufSize);
    toProcSize << strSize;

    UOPstream toProc(node, pBufGrid);
    toProc.write(ostr.str().data(), strSize);
} // sendGrids


template<typename GridType>
void Foam::foamVDB::sendGrid
(
    const typename GridType::Ptr grid,
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    std::vector<typename GridType::Ptr> inputGrids(1, grid);

    sendGrids<GridType>
    (
        inputGrids,
        pBufSize,
        pBufGrid,
        node
    );
}


template<typename GridType>
std::vector<typename GridType::Ptr> Foam::foamVDB::receiveGrids
(
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    // non-blocking communication
    UIPstream fromProcSize(node, pBufSize);
    std::streamsize strSize;
    fromProcSize >> strSize;

    char * incoming = new char[strSize];

    UIPstream fromProc(node, pBufGrid);
    fromProc.read(incoming, strSize);

    std::stringstream ss;
    ss.write(incoming, strSize);

    std::istringstream istr;
    istr.str(ss.str());

    openvdb::io::Stream strm(istr);
    openvdb::GridPtrVecPtr streamGrids = strm.getGrids();

    std::vector<typename GridType::Ptr> grids;

    for (openvdb::GridPtrVec::iterator it = streamGrids->begin(); it < streamGrids->end(); ++it)
    {
        typename GridType::Ptr grid = openvdb::gridPtrCast<GridType>(*it);

        grids.push_back(grid);
    }

    return grids;
} // receiveGrids


template<typename GridType>
typename GridType::Ptr Foam::foamVDB::receiveGrid
(
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    return receiveGrids<GridType>
    (
        pBufSize,
        pBufGrid,
        node
    )[0];
}


template<typename GridType>
void Foam::foamVDB::combineGrids
(
    const word name,
    std::vector<std::vector<typename GridType::Ptr>>& cellLevelGridsCollection,
    std::vector<typename GridType::Ptr>& cellLevelGrids
)
{
    Timer timer("combineGrids (topology union) " + name);

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; cellLevel++)
    {
        uLabel masterSet = 0;

        // find first valid set to use as master
        for (size_t setI = 0; setI < cellLevelGridsCollection.size(); setI++)
        {
            //check if valid pointer
            if (cellLevelGridsCollection[setI][cellLevel])
            {
                masterSet = setI;
                break;
            }
        }

        if (!cellLevelGridsCollection[masterSet][cellLevel])
        {
            cellLevelGridsCollection[masterSet][cellLevel] =
                GridType::create();
        }

        for (size_t setI = 0; setI < cellLevelGridsCollection.size(); setI++)
        {
            if (!cellLevelGridsCollection[setI][cellLevel]) continue;

            if (setI != masterSet)
            {
                cellLevelGridsCollection[masterSet][cellLevel]->topologyUnion(*cellLevelGridsCollection[setI][cellLevel]);

                // free memory, deletes managed object
                cellLevelGridsCollection[setI][cellLevel].reset();
            }
        }

        if (isValid<GridType>(cellLevelGrids[cellLevel]))
        {
            cellLevelGrids[cellLevel]->topologyUnion(*cellLevelGridsCollection[masterSet][cellLevel]);
        }
        else
        {
            cellLevelGrids[cellLevel] = cellLevelGridsCollection[masterSet][cellLevel];

            const word gridName = name + " level " + std::to_string(cellLevel);

            cellLevelGrids[cellLevel]->setName(gridName);
        }
    } // for cellLevel
} // combineGrids


template<class T>
struct SendListOp
{
    Foam::PstreamBuffers& pBuffers_;
    const Foam::List<T>& list_;

    SendListOp
    (
        Foam::PstreamBuffers& pBuffers,
        const Foam::List<T>& list
    )
    :
        pBuffers_(pBuffers),
        list_(list)
    {}

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t proci = range.begin(); proci < range.end(); ++proci)
        {
            size_t myProcNo = Foam::Pstream::myProcNo();

            if (proci == myProcNo) continue;

            Foam::UOPstream toProcI(proci, pBuffers_);

            toProcI << list_[myProcNo];
        }
    }
}; // SendListOp

template<class T>
struct ReceiveListOp
{
    Foam::PstreamBuffers& pBuffers_;
    Foam::List<T>& list_;

    ReceiveListOp
    (
        Foam::PstreamBuffers& pBuffers,
        Foam::List<T>& list
    )
    :
        pBuffers_(pBuffers),
        list_(list)
    {}

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        for (size_t proci = range.begin(); proci < range.end(); ++proci)
        {
            size_t myProcNo = Foam::Pstream::myProcNo();

            if (proci == myProcNo) continue;

            Foam::UIPstream fromProcI(proci, pBuffers_);

            T in(fromProcI);

            list_[proci] = in;
        }
    }
}; // ReceiveListOp

template<class T>
void Foam::foamVDB::allGatherList
(
    List<T>& list
)
{
    // threads are enabled by executable!
    bool threaded = Pstream::haveThreads();

    PstreamBuffers pBuffers(Pstream::commsTypes::nonBlocking);

    // send to all procs
    if (threaded)
    {
        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, Pstream::nProcs()),
            SendListOp<T>
            (
                pBuffers,
                list
            )
        );
    }
    else
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci == Pstream::myProcNo()) continue;

            UOPstream toProcI(proci, pBuffers);

            toProcI << list[Pstream::myProcNo()];
        }
    }

    pBuffers.finishedSends();

    // receive from all procs
    if (threaded)
    {
        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, Pstream::nProcs()),
            ReceiveListOp<T>
            (
                pBuffers,
                list
            )
        );
    }
    else
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci == Pstream::myProcNo()) continue;

            UIPstream fromProcI(proci, pBuffers);

            T in(fromProcI);

            list[proci] = in;
        }
    }
} // allGatherList


template<typename GridType>
Foam::List<Foam::List<Foam::boundBox>> Foam::foamVDB::sendReceiveBoundBoxes
(
    const std::vector<typename GridType::Ptr>& grids
)
{
    std::string gridName = grids[0]->getName();
    Timer timer("sendAndReceiveBoundBoxes " + gridName);

    // boundBox for each cellLevel
    List<boundBox> boundBoxes(grids.size());

    // populate boundBoxes for myProcNo
    for (size_t cellLevel = 0; cellLevel < grids.size(); ++cellLevel)
    {
        openvdb::CoordBBox bb = grids[cellLevel]->evalActiveVoxelBoundingBox();

        boundBoxes[cellLevel] =
            boundBox
            (
                point
                (
                    bb.min().x(),
                    bb.min().y(),
                    bb.min().z()
                ),
                point
                (
                    bb.max().x() + 1,
                    bb.max().y() + 1,
                    bb.max().z() + 1
                )
            );
    }

    List<List<boundBox>> procBoundBoxes(Pstream::nProcs());

    procBoundBoxes[Pstream::myProcNo()] = boundBoxes;

    // share with other processors
    allGatherList<List<boundBox>>(procBoundBoxes);

    return procBoundBoxes;
} //sendReceiveBoundBoxes


template<class Type>
void Foam::foamVDB::createVolField
(
    fvMesh& mesh,
    const std::vector<Type>& field,
    const word& name
)
{
    Timer timer("createVolField " + name);

    //Info<<"Creating field " << name << endl;

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    fieldType* volFieldPtr
    (
        new fieldType
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensioned<Type>(name, dimless, Zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );

    fieldType& volField = *volFieldPtr;

    //GGG
    bool threaded = true;
    if (threaded)
    {
        const auto grainSize =
            std::max<size_t>
            (
                volField.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        const tbb::blocked_range<size_t> range(0, volField.size(), grainSize);

        tbb::parallel_for
        (
            range,
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < range.end(); ++celli)
                {
                    volField[celli] = field[celli];
                }
            },
            tbb::simple_partitioner()
        );
    }
    else
    {
        forAll(volField, celli)
        {
            volField[celli] = field[celli];
        }
    }

    mesh.objectRegistry::store(volFieldPtr);
} // createVolField


template<typename GridType>
Foam::label Foam::foamVDB::getActiveVoxels
(
    const std::vector<typename GridType::Ptr>& grids
)
{
    Timer t("getActiveVoxels");

    label activeVoxels = 0;

    //TODO parallel_reduce + leafManager
    for (size_t i = 0; i < grids.size(); i++)
    {
        activeVoxels += grids[i]->activeVoxelCount();
    }

    return activeVoxels;
}

// ************************************************************************* //
