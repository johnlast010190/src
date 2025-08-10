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
    (c) 2021 ESI Ltd.

\*---------------------------------------------------------------------------*/

#include "parhipDecomp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "fields/Fields/Field/SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parhipDecomp, 0);

    addToRunTimeSelectionTable(decompositionMethod, parhipDecomp, dictionary);

    template<>
    const char* NamedEnum<parhipDecomp::configs, 6>::
    names[] =
    {
        "fast",
        "eco",
        "strong",
        "fast-social",
        "eco-social",
        "strong-social"
    };
}

const Foam::NamedEnum<Foam::parhipDecomp::configs, 6>
    Foam::parhipDecomp::configNames;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parhipDecomp::writeGraphSimple
(
    const fileName& meshPath,
    idxtype* vtxdist,
    idxtype* xadj,
    idxtype* adjncy,
    MPI_Comm* comm
) const
{
    word graphFile(meshPath.name()+".graph");

    std::ofstream str(graphFile.c_str());

    Info<< "Dumping KaHIP graph file to " << graphFile << endl
        << "Use this in combination with graph_checker and graph2binary." << endl;

    int rank, size;
    MPI_Comm_rank( *comm, &rank);
    MPI_Comm_size( *comm, &size);

    idxtype local_number_of_nodes = vtxdist[rank+1] - vtxdist[rank];
    idxtype local_number_of_edges = xadj[local_number_of_nodes];
    idxtype number_of_nodes = vtxdist[size];

    idxtype global_number_of_edges = 0;

    MPI_Allreduce(&local_number_of_edges, &global_number_of_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, *comm);

    if (Pstream::master())
    {
        str << std::to_string(number_of_nodes) << " " << std::to_string(global_number_of_edges/2) << std::endl;

        for (size_t celli = 0; celli < local_number_of_nodes; ++celli)
        {
            const label start = xadj[celli];
            const label end = xadj[celli+1];

            for (label i = start; i < end; ++i)
            {
                str << std::to_string(adjncy[i]+1) << " ";
            }
            str << "\n";
        }
        str.close();
    }

    for (label proci = Pstream::firstSlave(); proci < Pstream::nProcs(); ++proci)
    {
        MPI_Barrier(*comm);

        if (proci == Pstream::myProcNo())
        {
            std::ofstream str;
            str.open(graphFile.c_str(), std::ofstream::out | std::ofstream::app);
            for (size_t celli = 0; celli < local_number_of_nodes; ++celli)
            {
                const label start = xadj[celli];
                const label end = xadj[celli+1];

                for (label i = start; i < end; ++i)
                {
                    str << std::to_string(adjncy[i]+1) << " ";
                }
                str << "\n";
            }
            str.close();
        }
    }
} //writeGraphSimple


// Call parhip with options from dictionary.
Foam::label Foam::parhipDecomp::decompose
(
    const fileName& meshPath,
    const List<label>& adjncy,
    const List<label>& xadj,
    const scalarField& cWeights,
    std::vector<idxtype>& finalDecomp
) const
{
    List<label> dummyAdjncy(1);
    List<label> dummyXadj(1);
    dummyXadj[0] = 0;

    return decompose
    (
        meshPath,
        adjncy.size(),
        (adjncy.size() ? adjncy.begin() : dummyAdjncy.begin()),
        xadj.size(),
        (xadj.size() ? xadj.begin() : dummyXadj.begin()),
        cWeights,
        finalDecomp
    );
}


// Call parhip with options from dictionary.
Foam::label Foam::parhipDecomp::decompose
(
    const fileName& meshPath,
    const label adjncySize,
    const label adjncy[],
    const label xadjSize,
    const label xadj[],
    const scalarField& cWeights,

    std::vector<idxtype>& finalDecomp
) const
{
    if (debug)
    {
        Pout<< "parhipDecomp : entering with xadj:" << xadjSize << endl;
    }

    // Default setup
    enum configs kahipConfig = configs::FAST;
    double imbalance = 0.01;
    int seed = 0;
    bool verbose = false;
    bool writeGraph = false;

    if (decompositionDict_.found("kahipCoeffs"))
    {
        const dictionary& kahipCoeffs =
            decompositionDict_.subDict("kahipCoeffs");

        if (kahipCoeffs.found("config"))
        {
            kahipConfig = configNames.read(kahipCoeffs.lookup("config"));
        }
        kahipCoeffs.readIfPresent("imbalance", imbalance);
        kahipCoeffs.readIfPresent("verbose", verbose);

        //List<label> processorWeights;
        //kahipCoeffs.readIfPresent("processorWeights", processorWeights);
        //
        //if (processorWeights.size())
        //{
        //    if (debug)
        //    {
        //        Info<< "parhipDecomp : Using processor weights "
        //            << processorWeights
        //            << endl;
        //    }
        //}
        //
        kahipCoeffs.readIfPresent("seed", seed);

        kahipCoeffs.readIfPresent("writeGraph", writeGraph);
    }

    Info<< "parhipDecomp :"
        << " config=" << configNames[kahipConfig]
        << " imbalance=" << imbalance
        << " seed=" << seed
        << nl << endl;

    label numCells = xadjSize-1;

    // Cell weights (so on the vertices of the dual)
    std::vector<idxtype> cellWeights;

    // Check for externally provided cellweights and if so initialise weights
    const scalar minWeights = gMin(cWeights);
    const scalar maxWeights = gMax(cWeights);

    if (maxWeights > minWeights)
    {
        if (minWeights <= 0)
        {
            WarningInFunction
                << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
                << exit(FatalError);
        }

        if (cWeights.size())
        {
            // Convert to idxtype.
            cellWeights.reserve(cWeights.size());
            for (label i = 0; i < cWeights.size(); i++)
            {
                cellWeights.emplace_back
                (
                    idxtype(cWeights[i]/minWeights)
                );
            }
        }
        else
        {
            // Locally zero cells but not globally. Make sure we have
            // some size so .begin() does not return null pointer. Data
            // itself is never used.
            cellWeights.resize(1);
            cellWeights[0] = 1;
        }
    }


    // Number of partitions
    int nParts = nProcessors_;

    // Output: number of cut edges
    int edgeCut = 0;

    // Output: cell -> processor addressing
    finalDecomp.reserve(numCells);

    labelList procCells(Pstream::nProcs(), 0);
    procCells[Pstream::myProcNo()] = numCells;
    IPstream::allGatherList(procCells);

    label vtxdistSize = Pstream::nProcs() + 1;
    std::vector<idxtype> vtxdist;
    vtxdist.reserve(vtxdistSize);
    vtxdist.emplace_back(0);
    for (label i = 1; i < vtxdistSize; i++)
    {
        vtxdist.emplace_back(vtxdist[i-1] + procCells[i-1]);
    }

    std::vector<idxtype> xadj_param;
    xadj_param.reserve(xadjSize);
    for (label i = 0; i < xadjSize; i++)
    {
        xadj_param.emplace_back(xadj[i]);
    }

    std::vector<idxtype> adjncy_param;
    adjncy_param.reserve(adjncySize);
    for (label i = 0; i < adjncySize; i++)
    {
        adjncy_param.emplace_back(adjncy[i]);
    }

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    if (writeGraph)
    {
        writeGraphSimple
        (
            meshPath,
            vtxdist.data(),
            xadj_param.data(),
            adjncy_param.data(),
            &comm
        );
    }

    ParHIPPartitionKWay
    (
        vtxdist.data(),
        xadj_param.data(),
        adjncy_param.data(),
        (cWeights.size() ? cellWeights.data() : nullptr),
        nullptr,
        &nParts,
        &imbalance,
        !verbose,
        seed,
        int(kahipConfig),
        &edgeCut,
        finalDecomp.data(),
        &comm
    );

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parhipDecomp::parhipDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::parhipDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh.nCells()
            << exit(FatalError);
    }

    bool transparentAMIs =
        decompositionDict_.lookupOrDefault("transparentAMIs", false);

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli


    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells,
        transparentAMIs
    );

    // Decompose using default weights
    std::vector<idxtype> finalDecomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(points.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


Foam::labelList Foam::parhipDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    if (agglom.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        agglom,
        agglomPoints.size(),
        true,
        cellCells
    );

    // Decompose using weights
    std::vector<idxtype> finalDecomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        finalDecomp
    );

    // Rework back into decomposition for original mesh
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::parhipDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells(globalCellCells);

    // Decompose using weights
    std::vector<idxtype> finalDecomp;
    decompose
    (
        "parhip",
        cellCells.m(),
        cellCells.offsets(),
        cWeights,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(cellCentres.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //
