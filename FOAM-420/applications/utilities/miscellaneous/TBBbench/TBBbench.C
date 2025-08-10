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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/
#include <cassert>

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tick_count.h>

#include "cfdTools/general/include/fvCFD.H"
#include "global/argList/argList.H"
#include "memInfo/memInfo.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "db/IOobjectList/IOobjectList.H"
#include "include/parallelForAll.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

using namespace Foam;

// * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label checkProcessorFolders(const fileName& dir)
{
    label nProcs = 0;

    while
    (
        isDir
        (
            dir/"processor"
          + Foam::name(++nProcs)
        )
    )
    {}

    //wait for all procs to finish counting
    gMax(labelList(Pstream::nProcs(), 0));

    return nProcs;

    //// create dummy source folders
    //if (nProcs < Pstream::nProcs())
    //{
    //    instantList timeDirs;
    //    const word procFolder = "processor" + Foam::name(Pstream::myProcNo());
    //    if (Pstream::master())
    //    {
    //        const bool oldParRun = Pstream::parRun();
    //        Pstream::parRun() = false;
    //        timeDirs = Time::findTimes(dir/procFolder, "constant");
    //        Pstream::parRun() = oldParRun;
    //    }
    //    Pstream::scatter(timeDirs);
    //    forAll(timeDirs, i)
    //    {
    //        mkDir(dir/procFolder/timeDirs[i].name());
    //    }
    //}
    //else if (nProcs > Pstream::nProcs())
    //{
    //    FatalError
    //        << "Running with " << Pstream::nProcs()
    //        << " processors but case "
    //        << dir << " is decomposed in "
    //        << nProcs << " processors. "
    //        << nl << "Please run with maximum number of source and target processors"
    //        << exit(FatalError);
    //}
}


template<class Type>
void appendFieldNames
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    wordList& names
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            Info<<"Appending " << fieldName << endl;

            names.append(fieldName);
        }
    }
}


inline void printTime(const word name, const scalar timeSerial, const scalar timeParallel = -1.0)
{
    Info<< setw(40) << name << " - elapsed time: "
        << setw(5)
        << (timeParallel > 0 ? timeParallel : timeSerial)
        << setprecision(3) << " seconds";

    if (timeParallel > 0)
    {
        const scalar speedup = timeSerial/timeParallel;

        const size_t nThreads =
            tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        Info<< " - Speedup: " << setw(4) << speedup << setprecision(2)
            << " - Efficiency: " << setw(4) << 100*speedup/nThreads << setprecision(3)
            << "% (on " << nThreads << " threads)"
            << endl;
    }

    Info<<endl;
}

inline void parallel_for_mutex_count_aff
(
    const labelList& list,
    labelList& ncfm,
    std::vector<tbb::spin_mutex>& mut,
    const label grainSize,
    tbb::affinity_partitioner& affinity
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                tbb::spin_mutex::scoped_lock  lock{mut[celli]};
                ncfm[celli]++;
            }
        },
        affinity
    );
}

template<typename PartitionerType>
inline void parallel_for_mutex_count
(
    const labelList& list,
    labelList& ncfm,
    std::vector<tbb::spin_mutex>& mut,
    const label grainSize
)
{
    PartitionerType p;

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                tbb::spin_mutex::scoped_lock  lock{mut[celli]};
                ncfm[celli]++;
            }
        },
        p
    );
}

inline void parallel_for_mutex_fill_aff
(
    const labelList& list,
    labelList& ncfm,
    cellList& cellFaceAddrMutex,
    std::vector<tbb::spin_mutex>& mut,
    const label grainSize,
    tbb::affinity_partitioner& affinity
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                tbb::spin_mutex::scoped_lock  lock{mut[celli]};
                cellFaceAddrMutex[celli][ncfm[celli]++] = facei;
            }
        },
        affinity
    );
}

template<typename PartitionerType>
inline void parallel_for_mutex_fill
(
    const labelList& list,
    labelList& ncfm,
    cellList& cellFaceAddrMutex,
    std::vector<tbb::spin_mutex>& mut,
    const label grainSize
)
{
    PartitionerType p;

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                tbb::spin_mutex::scoped_lock  lock{mut[celli]};
                cellFaceAddrMutex[celli][ncfm[celli]++] = facei;
            }
        },
        p
    );
}

inline void parallel_for_atomic_count
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size()),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                ncfa[celli]++;
            }
        }
    );
}

inline void parallel_for_atomic_count
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa,
    const label grainSize
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                ncfa[celli]++;
            }
        },
        tbb::simple_partitioner()
    );
}

inline void parallel_for_atomic_count_aff
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa,
    const label grainSize,
    tbb::affinity_partitioner& affinity
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                ncfa[celli]++;
            }
        },
        affinity
    );
}

inline void parallel_for_atomic_fill
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa,
    cellList& cellFaceAddrAtomic
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size()),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                cellFaceAddrAtomic[celli][ncfa[celli]++] = facei;
            }
        }
    );
}

inline void parallel_for_atomic_fill
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa,
    cellList& cellFaceAddrAtomic,
    const label grainSize
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                cellFaceAddrAtomic[celli][ncfa[celli]++] = facei;
            }
        },
        tbb::simple_partitioner()
    );
}

inline void parallel_for_atomic_fill_aff
(
    const labelList& list,
    std::vector<std::atomic<label>>& ncfa,
    cellList& cellFaceAddrAtomic,
    const label grainSize,
    tbb::affinity_partitioner& affinity
)
{
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, list.size(), grainSize),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t facei = r.begin(); facei < r.end(); facei++)
            {
                const label& celli = list[facei];
                cellFaceAddrAtomic[celli][ncfa[celli]++] = facei;
            }
        },
        affinity
    );
}


// very lightweight struct for reading basic mesh data
// and multi-threaded calculation of face-cell addressing
struct polyMeshBasic
{
    pointIOField      points_;
    faceCompactIOList faces_;
    labelIOList       owner_;
    labelIOList       neighbour_;
    boundBox          bounds_;
    label             nCells_;
    cellList          cellFaceAddr_;
    label             grainSize_;

    polyMeshBasic
    (
        const Time& runTime,
        const fileName& pointsInstance,
        const fileName& facesInstance,
        const label grainSize = 1
    )
    :
        points_
        (
            IOobject
            (
                "points",
                pointsInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        faces_
        (
            IOobject
            (
                "faces",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        owner_
        (
            IOobject
            (
                "owner",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        neighbour_
        (
            IOobject
            (
                "neighbour",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        bounds_(points_, /*doReduce*/false),
        nCells_(-1),
        cellFaceAddr_(),
        grainSize_(grainSize)
    {}

    inline label nCells() const
    {
        return nCells_;
    }

    const boundBox& bounds() const
    {
        return bounds_;
    }

    const pointField& points() const
    {
        return points_;
    }

    const faceList& faces() const
    {
        return faces_;
    }

    const cellList& cells() const
    {
        return cellFaceAddr_;
    }

    struct ReduceCalcCells
    {
        const labelList& llist_;
        Map<labelList> faceCell_;

        ReduceCalcCells
        (
            const labelList& llist
        )
        :
            llist_(llist),
            faceCell_()
        {}

        ReduceCalcCells
        (
            const ReduceCalcCells& rhs,
            tbb::split
        )
        :
            llist_(rhs.llist_),
            faceCell_()
        {}

        void operator()(const tbb::blocked_range<size_t>& r)
        {
            //if (Pstream::master())
            //{
            //    std::cout<< "    reduce thread ID "
            //        << std::this_thread::get_id()
            //        << " - range " << r.begin()
            //        << ":" << r.end()
            //        << std::endl;
            //}
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                const label& celli = llist_[i];
                faceCell_(celli).append(i);
            }
        }

        void join(const ReduceCalcCells& rhs)
        {
            auto keyIterPair = rhs.faceCell_.keys();
            for (const auto& i : keyIterPair)
            {
                const labelList& rhsFaces = rhs.faceCell_[i];
                faceCell_(i).append(rhsFaces);
            }
        }
    }; //ReduceCalcCells


    struct ReduceBody
    {
        const labelList& llist_;
        label nCells_;

        ReduceBody
        (
            const labelList& llist,
            const label init
        )
        :
            llist_(llist),
            nCells_(init)
        {}

        ReduceBody
        (
            const ReduceBody& rhs,
            tbb::split
        )
        :
            llist_(rhs.llist_),
            nCells_(rhs.nCells_)
        {}

        void operator()(const tbb::blocked_range<size_t>& r)
        {
            //if (Pstream::master())
            //{
            //    std::cout<< "    reduce thread ID "
            //        << std::this_thread::get_id()
            //        << " - range " << r.begin()
            //        << ":" << r.end()
            //        << std::endl;
            //}
            // Performance. Sometimes putting frequently accessed
            // values into local variables helps the compiler
            // optimize the loop better, because local variables
            // are often easier for the compiler to track.
            label n = nCells_;

            const labelList& myList = llist_;

            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                n = std::max(n, myList[i]);
            }

            nCells_ = n;
        }

        void join(const ReduceBody& rhs)
        {
            nCells_ = std::max(nCells_, rhs.nCells_);
        }
    }; //ReduceBody

    label reduceBody(const labelList& llist, label identity, label grainSize)
    {
        ReduceBody func(llist, identity);

        const tbb::blocked_range<size_t> range(0, llist.size(), grainSize);

        tbb::parallel_reduce(range, func, tbb::simple_partitioner());

        return func.nCells_;
    }

    label reduce(const labelList& llist, label identity, label grainSize)
    {
        const tbb::blocked_range<size_t> range(0, llist.size(), grainSize);

        return tbb::parallel_reduce
        (
            range,
            identity,
            ///*func*/[&](const tbb::blocked_range<size_t>& r, label init) -> label
            /*func*/[&](const tbb::blocked_range<size_t>& r, const label& init)
            {
                //if (Pstream::master())
                //{
                //    std::cout<< "    reduce thread ID "
                //        << std::this_thread::get_id()
                //        << " - range " << r.begin()
                //        << ":" << r.end()
                //        << std::endl;
                //}
                label res = init;
                for (size_t i = r.begin(); i != r.end(); ++i)
                {
                    //init = Foam::max(init, llist[i]);
                    res = Foam::max(res, llist[i]);
                }
                //return init;
                return res;
            },
            ///*reduction*/[](label x, label y) -> label
            /*reduction*/[](const label& x, const label& y)
            {
                return Foam::max(x, y);
            },
            //tbb::static_partitioner()
            tbb::simple_partitioner()
        );
    }

    label countCells(const labelList& own, const labelList& nei, const label grainSize)
    {
        //SERIAL
        tbb::tick_count t0S = tbb::tick_count::now();
        label nCellsS = -1;

        forAll(own, facei)
        {
            nCellsS = Foam::max(nCellsS, own[facei]);
        }
        forAll(nei, facei)
        {
            nCellsS = Foam::max(nCellsS, nei[facei]);
        }
        ++nCellsS;

        scalar tSerial = (tbb::tick_count::now() - t0S).seconds();
        printTime("    countCells serial", tSerial);

        Info<< "    nCells (processor0) = " << nCellsS << endl;
        Info<< "    own.size() " << own.size() << " - nei.size() " << nei.size() << endl;
        Info<< "    avg time/iter " << tSerial/((own.size()+nei.size())/2.0) << endl;

        Info<<endl;

        {
        //////////////
        // TASKS
        tbb::tick_count t0T = tbb::tick_count::now();
        tbb::task_group tasks;
        label maxOwn = -1, maxNei = -1;
        tasks.run
        (
            [&maxOwn, &own]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                forAll(own, facei)
                {
                    maxOwn = Foam::max(maxOwn, own[facei]);
                }
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&maxNei, &nei]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                forAll(nei, facei)
                {
                    maxNei = Foam::max(maxNei, nei[facei]);
                }
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2", (t1-t0).seconds());
            }
        );
        tasks.wait();
        label nCellsT = Foam::max(maxNei, maxOwn) + 1;

        assert(nCellsT == nCellsS);

        scalar tTasks = (tbb::tick_count::now() - t0T).seconds();
        printTime("    countCells tasks", tSerial, tTasks);

        Info<<endl;
        }

        {
        //////////////
        // TASKS run_and_wait
        tbb::tick_count t0T = tbb::tick_count::now();
        tbb::task_group tasks;
        label maxOwn = -1, maxNei = -1;
        tasks.run
        (
            [&maxOwn, &own]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                forAll(own, facei)
                {
                    maxOwn = Foam::max(maxOwn, own[facei]);
                }
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1", (t1-t0).seconds());
            }
        );
        tasks.run_and_wait
        (
            [&maxNei, &nei]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                forAll(nei, facei)
                {
                    maxNei = Foam::max(maxNei, nei[facei]);
                }
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2", (t1-t0).seconds());
            }
        );
        label nCellsT = Foam::max(maxNei, maxOwn) + 1;

        assert(nCellsT == nCellsS);

        scalar tTasks = (tbb::tick_count::now() - t0T).seconds();
        printTime("    countCells tasks run_and_wait", tSerial, tTasks);

        Info<<endl;
        }

        //////////////
        //TBB parallel_reduce lambda
        {
        tbb::tick_count t0L = tbb::tick_count::now();

        label nCellsL = reduce(own, -1, grainSize);

        tbb::tick_count t1L = tbb::tick_count::now();
        printTime("        lambda reduce1", (t1L-t0L).seconds());

        nCellsL = reduce(nei, nCellsL, grainSize);
        ++nCellsL;

        tbb::tick_count t2L = tbb::tick_count::now();
        printTime("        lambda reduce2", (t2L-t1L).seconds());

        scalar tLambda = (t2L - t0L).seconds();
        printTime("    countCells reduce(lambda)", tSerial, tLambda);

        assert(nCellsL == nCellsS);
        Info<<endl;
        }

        //////////////
        //TBB parallel_reduce with Body
        {
        tbb::tick_count t0B = tbb::tick_count::now();
        label nCellsB = reduceBody(own, -1, grainSize);

        tbb::tick_count t1B = tbb::tick_count::now();
        printTime("        body reduce1", (t1B-t0B).seconds());

        nCellsB = reduceBody(nei, nCellsB, grainSize);
        ++nCellsB;

        tbb::tick_count t2B = tbb::tick_count::now();
        printTime("        body reduce2", (t2B-t1B).seconds());

        scalar tBody = (t2B - t0B).seconds();
        printTime("    countCells reduce(body)", tSerial, tBody);

        assert(nCellsB == nCellsS);
        }

        return nCellsS;
    }

    void calcCells(const label grainSize)
    {
        Info<< "\n    ======== COUNT CELLS =======" << nl <<endl;
        nCells_ = countCells(owner_, neighbour_, grainSize);

        Info<< "\n    ======= CALC ADDRESSING =======" << nl <<endl;
        calcCells(cellFaceAddr_, owner_, neighbour_, nCells_, grainSize);

        Info<< endl;
    }

    class FaceCellAddr
    {
        cellList& cellFaceAddr_;
        const labelList& own_;
        const labelList& nei_;

    public:

        FaceCellAddr
        (
            cellList& cellFaceAddr,
            const labelList& own,
            const labelList& nei
        )
        :
            cellFaceAddr_(cellFaceAddr),
            own_(own),
            nei_(nei)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {

            labelList ncf(range.size(), 0);

            label ownStart = -1, ownEnd = -1;
            label neiStart = -1, neiEnd = -1;

            forAll(own_, facei)
            {
                const label& celli = own_[facei];

                if (celli >= range.begin() && celli < range.end())
                {
                    ownEnd = facei;
                    if (ownStart == -1) ownStart = facei;

                    const label idx = celli - range.begin();

                    ncf[idx]++;
                }
            }

            forAll(nei_, facei)
            {
                const label& celli = nei_[facei];

                if (celli >= range.begin() && celli < range.end())
                {
                    neiEnd = facei;
                    if (neiStart == -1) neiStart = facei;

                    const label idx = celli - range.begin();

                    ncf[idx]++;
                }
            }
            //if (Pstream::master())
            //{
            //    std::cout<< "    FaceCell thread ID "
            //        << std::this_thread::get_id()
            //        << " - range " << range.begin()
            //        << ":" << range.end()
            //        << " - size " << range.size()
            //        << " - ownRange " << ownStart << ":" << ownEnd
            //        << " - neiRange " << neiStart << ":" << neiEnd
            //        << std::endl;
            //}

            forAll(ncf, i)
            {
                const label celli = i + range.begin();

                cellFaceAddr_[celli].setSize(ncf[i]);
                ncf[i] = 0;
            }
            //tbb::parallel_for
            //(
            //    tbb::blocked_range<size_t>(0, ncf.size()),
            //    [&](const tbb::blocked_range<size_t>& r)
            //    {
            //        for (size_t i = r.begin(); i < r.end(); i++)
            //        {
            //            const label celli = i + range.begin();
            //
            //            cellFaceAddr_[celli].setSize(ncf[i]);
            //            ncf[i] = 0;
            //        }
            //    }
            //);

            //forAll(own_, facei)
            for (label facei = ownStart; facei <= ownEnd; facei++)
            {
                const label& celli = own_[facei];

                if (celli >= range.begin() && celli < range.end())
                {
                    const label idx = celli - range.begin();

                    cellFaceAddr_[celli][ncf[idx]++] = facei;
                }
            }

            //forAll(nei_, facei)
            for (label facei = neiStart; facei <= neiEnd; facei++)
            {
                const label& celli = nei_[facei];

                if (celli >= range.begin() && celli < range.end())
                {
                    const label idx = celli - range.begin();

                    cellFaceAddr_[celli][ncf[idx]++] = facei;
                }
            }
        }
    }; //class FaceCellAddr

    void calcCells
    (
        cellList& cellFaceAddr,
        const labelList& own,
        const labelList& nei,
        const label nCells,
        const label grainSize
    )
    {
        //SERIAL
        tbb::tick_count t0S = tbb::tick_count::now();
        cellList cellFaceAddrSerial;
        cellFaceAddrSerial.setSize(nCells);

        labelList ncf(nCells, 0);

        forAll(own, facei)
        {
            ncf[own[facei]]++;
        }

        forAll(nei, facei)
        {
            if (nei[facei] >= 0)
            {
                ncf[nei[facei]]++;
            }
        }

        forAll(cellFaceAddrSerial, celli)
        {
            cellFaceAddrSerial[celli].setSize(ncf[celli]);
        }
        ncf = 0;

        forAll(own, facei)
        {
            label celli = own[facei];

            cellFaceAddrSerial[celli][ncf[celli]++] = facei;
        }

        forAll(nei, facei)
        {
            label celli = nei[facei];

            if (celli >= 0)
            {
                cellFaceAddrSerial[celli][ncf[celli]++] = facei;
            }
        }
        scalar tSerial = (tbb::tick_count::now() - t0S).seconds();
        printTime("    calcCells serial", tSerial);
        Info<< "    avg time/iter " << tSerial/((own.size()+nei.size())/2.0) << endl;
        Info<<endl;

        //{
        ////////////////
        //// parallel_for
        //tbb::tick_count t0P = tbb::tick_count::now();

        //cellList cellFaceAddrPar;
        //cellFaceAddrPar.setSize(nCells);

        //const tbb::blocked_range<size_t> cellRange(0, nCells, grainSize);

        //tbb::parallel_for
        //(
        //    cellRange,
        //    FaceCellAddr(cellFaceAddrPar, own, nei),
        //    tbb::simple_partitioner()
        //);

        //scalar tPar = (tbb::tick_count::now() - t0P).seconds();
        //printTime("    calcCells parallel_for", tSerial, tPar);

        //assert(cellFaceAddrPar == cellFaceAddrSerial);
        //}

        //{
        ////////////////
        //// parallel_reduce
        ////////////////
        //tbb::tick_count t0R = tbb::tick_count::now();

        //cellList cellFaceAddrRed(nCells);

        //ReduceCalcCells ownCells(own);
        //tbb::parallel_reduce
        //(
        //    tbb::blocked_range<size_t>(0, own.size(), grainSize),
        //    ownCells,
        //    tbb::simple_partitioner()
        //);

        //scalar tRed = (tbb::tick_count::now() - t0R).seconds();
        //printTime("    calcCells parallel_reduce (partial)", tRed);
        ////assert(cellFaceAddrRed == cellFaceAddrSerial);
        //Info<<endl;
        //}

        {
        //////////////
        // spin_mutex simple_partitioner
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        parallel_for_mutex_count<tbb::simple_partitioner>(own, ncfm, mut, grainSize);
        parallel_for_mutex_count<tbb::simple_partitioner>(nei, ncfm, mut, grainSize);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                    ncfm[celli] = 0;
                }
            },
            tbb::simple_partitioner()
        );

        parallel_for_mutex_fill<tbb::simple_partitioner>(own, ncfm, cellFaceAddrMutex, mut, grainSize);
        parallel_for_mutex_fill<tbb::simple_partitioner>(nei, ncfm, cellFaceAddrMutex, mut, grainSize);

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells spin_mutex simple_partitioner", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // spin_mutex auto_partitioner
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        parallel_for_mutex_count<tbb::auto_partitioner>(own, ncfm, mut, grainSize);
        parallel_for_mutex_count<tbb::auto_partitioner>(nei, ncfm, mut, grainSize);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                    ncfm[celli] = 0;
                }
            },
            tbb::auto_partitioner()
        );

        parallel_for_mutex_fill<tbb::auto_partitioner>(own, ncfm, cellFaceAddrMutex, mut, grainSize);
        parallel_for_mutex_fill<tbb::auto_partitioner>(nei, ncfm, cellFaceAddrMutex, mut, grainSize);

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells spin_mutex auto_partitioner", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // spin_mutex static_partitioner
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        parallel_for_mutex_count<tbb::static_partitioner>(own, ncfm, mut, grainSize);
        parallel_for_mutex_count<tbb::static_partitioner>(nei, ncfm, mut, grainSize);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                    ncfm[celli] = 0;
                }
            },
            tbb::static_partitioner()
        );

        parallel_for_mutex_fill<tbb::static_partitioner>(own, ncfm, cellFaceAddrMutex, mut, grainSize);
        parallel_for_mutex_fill<tbb::static_partitioner>(nei, ncfm, cellFaceAddrMutex, mut, grainSize);

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells spin_mutex static_partitioner", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // spin_mutex affinity_partitioner OwnNei
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        static tbb::affinity_partitioner affOwn;
        static tbb::affinity_partitioner affNei;

        parallel_for_mutex_count_aff(own, ncfm, mut, grainSize, affOwn);
        parallel_for_mutex_count_aff(nei, ncfm, mut, grainSize, affNei);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                    ncfm[celli] = 0;
                }
            },
            affOwn
        );

        parallel_for_mutex_fill_aff(own, ncfm, cellFaceAddrMutex, mut, grainSize, affOwn);
        parallel_for_mutex_fill_aff(nei, ncfm, cellFaceAddrMutex, mut, grainSize, affNei);

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells spin_mutex affinity_partitioner OwnNei", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // spin_mutex affinity_partitioner
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        static tbb::affinity_partitioner aff;

        parallel_for_mutex_count_aff(own, ncfm, mut, grainSize, aff);
        parallel_for_mutex_count_aff(nei, ncfm, mut, grainSize, aff);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                    ncfm[celli] = 0;
                }
            },
            aff
        );

        parallel_for_mutex_fill_aff(own, ncfm, cellFaceAddrMutex, mut, grainSize, aff);
        parallel_for_mutex_fill_aff(nei, ncfm, cellFaceAddrMutex, mut, grainSize, aff);

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells spin_mutex affinity_partitioner unique", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // spin_mutex + tasks
        //////////////
        tbb::tick_count t0M = tbb::tick_count::now();

        cellList cellFaceAddrMutex;
        cellFaceAddrMutex.setSize(nCells);

        labelList ncfm(nCells, 0);

        std::vector<tbb::spin_mutex> mut(nCells);

        tbb::task_group tasks;
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_mutex_count<tbb::simple_partitioner>(own, ncfm, mut, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 count", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_mutex_count<tbb::simple_partitioner>(nei, ncfm, mut, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 count", (t1-t0).seconds());
            }
        );
        tasks.wait();

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrMutex[celli].setSize(ncfm[celli]);
                }
            },
            tbb::simple_partitioner()
        );
        ncfm = 0;

        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_mutex_fill<tbb::simple_partitioner>(own, ncfm, cellFaceAddrMutex, mut, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 fill", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_mutex_fill<tbb::simple_partitioner>(nei, ncfm, cellFaceAddrMutex, mut, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 fill", (t1-t0).seconds());
            }
        );
        tasks.wait();

        scalar tMut = (tbb::tick_count::now() - t0M).seconds();
        printTime("    calcCells nested (spin_mutex + tasks)", tSerial, tMut);
        assert(cellFaceAddrMutex == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // atomic
        //////////////
        tbb::tick_count t0A = tbb::tick_count::now();

        cellList cellFaceAddrAtomic;
        cellFaceAddrAtomic.setSize(nCells);

        std::vector<std::atomic<label>> ncfa(nCells);

        parallel_for_atomic_count(own, ncfa, grainSize);
        parallel_for_atomic_count(nei, ncfa, grainSize);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrAtomic[celli].setSize(ncfa[celli]);
                    ncfa[celli] = 0;
                }
            },
            tbb::simple_partitioner()
        );

        parallel_for_atomic_fill(own, ncfa, cellFaceAddrAtomic, grainSize);
        parallel_for_atomic_fill(nei, ncfa, cellFaceAddrAtomic, grainSize);

        scalar tAtom = (tbb::tick_count::now() - t0A).seconds();
        printTime("    calcCells atomic", tSerial, tAtom);
        assert(cellFaceAddrAtomic == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // atomic auto_partitioner
        //////////////
        tbb::tick_count t0A = tbb::tick_count::now();

        cellList cellFaceAddrAtomic;
        cellFaceAddrAtomic.setSize(nCells);

        std::vector<std::atomic<label>> ncfa(nCells);

        parallel_for_atomic_count(own, ncfa);
        parallel_for_atomic_count(nei, ncfa);

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrAtomic[celli].setSize(ncfa[celli]);
                    ncfa[celli] = 0;
                }
            }
        );

        parallel_for_atomic_fill(own, ncfa, cellFaceAddrAtomic);
        parallel_for_atomic_fill(nei, ncfa, cellFaceAddrAtomic);

        scalar tAtom = (tbb::tick_count::now() - t0A).seconds();
        printTime("    calcCells atomic auto_partitioner", tSerial, tAtom);
        assert(cellFaceAddrAtomic == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // atomic + tasks
        //////////////
        tbb::tick_count t0A = tbb::tick_count::now();

        cellList cellFaceAddrAtomic;
        cellFaceAddrAtomic.setSize(nCells);

        std::vector<std::atomic<label>> ncfa(nCells);

        tbb::task_group tasks;
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count(own, ncfa, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 count", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count(nei, ncfa, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 count", (t1-t0).seconds());
            }
        );
        tasks.wait();

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrAtomic[celli].setSize(ncfa[celli]);
                    ncfa[celli] = 0;
                }
            },
            tbb::simple_partitioner()
        );

        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill(own, ncfa, cellFaceAddrAtomic, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 fill", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill(nei, ncfa, cellFaceAddrAtomic, grainSize);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 fill", (t1-t0).seconds());
            }
        );
        tasks.wait();

        scalar tAtom = (tbb::tick_count::now() - t0A).seconds();
        printTime("    calcCells nested (atomic + tasks)", tSerial, tAtom);
        assert(cellFaceAddrAtomic == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // atomic + tasks auto_partitioner
        //////////////
        tbb::tick_count t0A = tbb::tick_count::now();

        cellList cellFaceAddrAtomic;
        cellFaceAddrAtomic.setSize(nCells);

        std::vector<std::atomic<label>> ncfa(nCells);

        tbb::task_group tasks;
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count(own, ncfa);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 count", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count(nei, ncfa);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 count", (t1-t0).seconds());
            }
        );
        tasks.wait();

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrAtomic[celli].setSize(ncfa[celli]);
                    ncfa[celli] = 0;
                }
            }
        );

        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill(own, ncfa, cellFaceAddrAtomic);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 fill", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill(nei, ncfa, cellFaceAddrAtomic);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 fill", (t1-t0).seconds());
            }
        );
        tasks.wait();

        scalar tAtom = (tbb::tick_count::now() - t0A).seconds();
        printTime("    calcCells nested (atomic + tasks) with auto_partitioner", tSerial, tAtom);
        assert(cellFaceAddrAtomic == cellFaceAddrSerial);

        Info<< endl;
        }

        {
        //////////////
        // atomic + tasks +affinity
        //////////////
        tbb::tick_count t0A = tbb::tick_count::now();

        cellList cellFaceAddrAtomic;
        cellFaceAddrAtomic.setSize(nCells);

        std::vector<std::atomic<label>> ncfa(nCells);

        static tbb::affinity_partitioner affOwn;
        static tbb::affinity_partitioner affNei;

        tbb::task_group tasks;
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count_aff(own, ncfa, grainSize, affOwn);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 count", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_count_aff(nei, ncfa, grainSize, affNei);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 count", (t1-t0).seconds());
            }
        );
        tasks.wait();

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells, grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddrAtomic[celli].setSize(ncfa[celli]);
                    ncfa[celli] = 0;
                }
            },
            affOwn
        );

        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill_aff(own, ncfa, cellFaceAddrAtomic, grainSize, affOwn);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task1 fill", (t1-t0).seconds());
            }
        );
        tasks.run
        (
            [&]
            {
                tbb::tick_count t0 = tbb::tick_count::now();
                parallel_for_atomic_fill_aff(nei, ncfa, cellFaceAddrAtomic, grainSize, affNei);
                tbb::tick_count t1 = tbb::tick_count::now();
                printTime("        task2 fill", (t1-t0).seconds());
            }
        );
        tasks.wait();

        scalar tAtom = (tbb::tick_count::now() - t0A).seconds();
        printTime("    calcCells nested (atomic + tasks) with affinity", tSerial, tAtom);
        assert(cellFaceAddrAtomic == cellFaceAddrSerial);

        Info<< endl;
        }

    } // calcCells
}; //polyMeshBasic


template<class Type>
void readFields
(
    const label fieldSize,
    const Time& runTime,
    const wordList& names,
    List<Field<Type>>& fields
)
{
    fields.setSize(names.size());

    // serial loop because read from stream
    // can be done by one thread at a time
    forAll(names, i)
    {
        const word& fieldName = names[i];

        localIOdictionary fieldDict
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            ),
            GeometricField<Type, fvPatchField, volMesh>::typeName
        );

        Field<Type> fieldSource("internalField", fieldDict, fieldSize);

        fields[i] = fieldSource;
    }
} // readFields


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "TBB benchmark"
    );

    argList::noFunctionObjects();

    argList::noCheckProcessorDirectories();

    // Create argList. This will check for non-existing processor dirs.
    // Create processor directory if non-existing
    Foam::argList args
    (
        argc,
        argv,
        /*checkArgs*/ true,
        /*checkOpts*/ true,
        /*initialise*/ true,
        /*needsThread*/true
    );

    Info<<"MPI_THREAD_MULTIPLE "
        << (Pstream::haveThreads() ? "true" : "false")
        << endl;

    //if (Pstream::parRun() && !isDir(args.path()))
    //{
    //    mkDir(args.path());
    //}

    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    Foam::Time runTimeTarget(Foam::Time::controlDictName, args);

    memInfo mem;
    Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    const word dictName("TBBbenchDict");

    Info<< "\nReading " << dictName << nl << endl;
    IOdictionary benchDict
    (
        IOobject
        (
            dictName,
            runTimeTarget.system(),
            "",
            runTimeTarget,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const labelList threads   = benchDict.lookupOrDefault<labelList>("threads", labelList(1, 0));
    const scalarList grainSize = benchDict.lookupOrDefault<scalarList>("grainSize", scalarList(1, 1.0));

    HashSet<word> selectedFields;
    if (benchDict.found("fields"))
    {
        benchDict.lookup("fields")() >> selectedFields;
    }

    const fileName rootDirSource = fileName(benchDict.lookup("sourceCase")).toAbsolute();

    if (!isDir(rootDirSource.path()))
    {
        FatalErrorInFunction
            << "source case directory: "
            << rootDirSource
            << " does not exist!"
            << exit(FatalError);
    }

    //check sourceCase processors
    const label nProcsSource = checkProcessorFolders(rootDirSource);

    Info<< "\nsourceCase " << rootDirSource
        << " decomposed in " << nProcsSource
        << " processors." << endl;

    word sourceTimeName("latestTime");

    if (benchDict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = benchDict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            sourceTimeName = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = readScalar(benchDict.lookup("sourceTime"));
            sourceTimeName = name(sourceTimeScalar);
        }
    }

    // assign a list of source processor folders to each node
    labelList nSourceMeshesInProc(Pstream::nProcs());
    labelList sourceMeshToProc(identity(nProcsSource));
    labelListList sourceMeshNoInProc =
        fvMeshDistribute::calcMeshToProcMap
        (
            nProcsSource,
            nSourceMeshesInProc,
            sourceMeshToProc
        );
    //e.g. 2 procs, 5 meshes
    // sourceMeshNoInProc
    // 2
    // (
    //   ( 0 2 3 )
    //   ( 1 4 )
    // )

    const labelList& mySourceMeshes =
        sourceMeshNoInProc[Pstream::myProcNo()];

    runTimeTarget.time().clockTimeIncrement();

    // initialize scalar and vector grids
    Time runTimeSource0
    (
        Time::controlDictName,
        rootDirSource,
        (
            nProcsSource > 1
          ? fileName(word("processor0"))
          : fileName(word("./"))
        )
    );

    instantList sourceTimes = runTimeSource0.times();
    label sourceTimeIndex = runTimeSource0.timeIndex();

    if (sourceTimeName == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(sourceTimeName);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTimeSource0.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    // Search for list of objects
    IOobjectList objects0(runTimeSource0, runTimeSource0.timeName());

    //Info<< "Objects found in sourceTime "
    //    << runTimeSource0.timeName()
    //    << objects0.toc() << endl;

    Info<< "\nReading source time "
        << runTimeSource0.timeName()
        << nl << endl;

    wordList scalarNames, vectorNames;

    appendFieldNames<scalar>
    (
         objects0,
         selectedFields,
         scalarNames
    );

    appendFieldNames<vector>
    (
         objects0,
         selectedFields,
         vectorNames
    );

    const fileName pointsInstance = runTimeSource0.findInstance(polyMesh::meshSubDir, "points");
    const fileName facesInstance  = runTimeSource0.findInstance(polyMesh::meshSubDir, "faces");

    Info<<"\nInitialized in "
        << runTimeTarget.time().clockTimeIncrement() << " (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << ") s" << nl << endl;

    // mutex for ISstream
    // only one thread at a time can read from stream
    tbb::spin_mutex readMutex;

    auto readMesh =
        [&](tbb::blocked_range<size_t>& range)
        {
            for (size_t i = range.begin(); i < range.end(); i++)
            {
                label proci = sourceMeshNoInProc[Pstream::myProcNo()][i];

                //if (Pstream::master())
                //{
                //    std::cout<< "readMesh thread ID "
                //        << tbb::this_tbb_thread::get_id()
                //        << std::endl;
                //}

                // lock ISstream for reading from disk
                tbb::spin_mutex::scoped_lock  streamLock;

                streamLock.acquire(readMutex);

                Time runTimeSourceI
                (
                    runTimeSource0.controlDict(),
                    rootDirSource,
                    (
                        nProcsSource > 1
                      ? fileName(word("processor") + name(proci))
                      : fileName(word("./"))
                    )
                );

                // assuming all processor folders have same times
                runTimeSourceI.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

                label nCells = 0;

                polyMeshBasic meshI(runTimeSourceI, pointsInstance, facesInstance, /*grainSizes*/1);

                Info<< "\nRead source mesh "
                    << (
                           nProcsSource > 1
                         ? rootDirSource/fileName(word("processor") + name(proci))
                         : rootDirSource
                       )
                    << " in "
                    << runTimeSourceI.time().clockTimeIncrement() << " (cpu "
                    << runTimeSourceI.time().cpuTimeIncrement() << ") s"
                    << nl << endl;

                // release lock
                streamLock.release();

                forAll(grainSize, gsI)
                {
                    const label gs = grainSize[gsI];

                    Info<< "  ~~~~~ Using grainSize " << gs << " ~~~~~" << endl;
                    meshI.calcCells(gs);

                    nCells = meshI.nCells();
                }

                // Read volFields
                // ~~~~~~~~~~~~~
                List<scalarField> scalarFields;
                List<vectorField> vectorFields;

                // do not share ISstream with other threads
                streamLock.acquire(readMutex);

                runTimeSourceI.time().clockTimeIncrement();
                runTimeSourceI.time().cpuTimeIncrement();

                readFields<scalar>
                (
                    nCells,
                    runTimeSourceI,
                    scalarNames,
                    scalarFields
                );

                readFields<vector>
                (
                    nCells,
                    runTimeSourceI,
                    vectorNames,
                    vectorFields
                );
                Info<< "\nRead fields of "
                    << (word("processor") + name(proci))
                    << " in "
                    << runTimeSourceI.time().clockTimeIncrement() << " (cpu "
                    << runTimeSourceI.time().cpuTimeIncrement() << ") s"
                    << nl << endl;

                streamLock.release();
            } //for i in range mySourceProc
        }; // readMesh


    // do not update processor boundaries
    bool oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    std::unique_ptr<tbb::global_control> control;

    forAll(threads, threadI)
    {
        const label nThreads = threads[threadI];

        if (nThreads > 0)
        {
            // note, threads == 0 means use all threads (default), so don't
            // manually create a tbb::global_control in this case
            control.reset
            (
                new tbb::global_control
                (
                    tbb::global_control::max_allowed_parallelism,
                    nThreads
                )
            );
        }

        const size_t nt = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        Info<< "\n~~~ Using " << nt << " thread"
            << (nt > 1 ? "s" : "")
            << " ~~~" << endl;

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>
            (
                0,
                mySourceMeshes.size()
            ),
            readMesh
        );
    } // for nThreads

    UPstream::parRun() = oldParRun;

    Info<< "\nelapsedCpuTime: " << runTimeTarget.elapsedCpuTime() << "s\n"
        << "elapsedClockTime: " << runTimeTarget.elapsedClockTime() << "s,\n"
        << mem.update().peak() << " kB (peak)"
        << nl << endl;

    return 0;
}

// ************************************************************************* //
