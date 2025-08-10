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
    (c) 2017-2021 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "metisLikeDecomp.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisLikeDecomp::decomposeGeneral
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    List<int>& decomp
) const
{
    if (!Pstream::parRun())
    {
        return decomposeSerial
        (
            adjncy,
            xadj,
            cWeights,
            decomp
        );
    }

    if (debug)
    {
        Info<< type() << "Decomp : running in parallel."
            << " Decomposing all of graph on master processor." << endl;
    }
    globalIndex globalCells(xadj.size()-1);
    label nTotalConnections = returnReduce(adjncy.size(), sumOp<label>());

    // Send all to master. Use scheduled to save some storage.
    if (Pstream::master())
    {
        List<label> allAdjncy(nTotalConnections);
        List<label> allXadj(globalCells.size()+1);
        List<scalar> allWeights(globalCells.size());

        // Insert my own
        label nTotalCells = 0;
        forAll(cWeights, celli)
        {
            allXadj[nTotalCells] = xadj[celli];
            allWeights[nTotalCells++] = cWeights[celli];
        }
        nTotalConnections = 0;
        forAll(adjncy, i)
        {
            allAdjncy[nTotalConnections++] = adjncy[i];
        }

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            List<label> nbrAdjncy(fromSlave);
            List<label> nbrXadj(fromSlave);
            List<scalar> nbrWeights(fromSlave);

            // Append.
            forAll(nbrXadj, celli)
            {
                allXadj[nTotalCells] = nTotalConnections+nbrXadj[celli];
                allWeights[nTotalCells++] = nbrWeights[celli];
            }
            // No need to renumber xadj since already global.
            forAll(nbrAdjncy, i)
            {
                allAdjncy[nTotalConnections++] = nbrAdjncy[i];
            }
        }
        allXadj[nTotalCells] = nTotalConnections;

        List<int> allDecomp;
        decomposeSerial
        (
            allAdjncy,
            allXadj,
            allWeights,
            allDecomp
        );


        // Send allFinalDecomp back
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            OPstream toSlave(Pstream::commsTypes::scheduled, slave);
            toSlave << SubList<int>
            (
                allDecomp,
                globalCells.localSize(slave),
                globalCells.offset(slave)
            );
        }

        // Get my own part (always first)
        decomp = SubList<int>(allDecomp, globalCells.localSize());
    }
    else
    {
        // Send my part of the graph (already in global numbering)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            toMaster
                << adjncy
                << SubList<label>(xadj, xadj.size()-1)
                << cWeights;
        }

        // Receive back decomposition
        IPstream fromMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );
        fromMaster >> decomp;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisLikeDecomp::metisLikeDecomp
(
    const word& derivedType,
    const dictionary& decompDict
)
:
    decompositionMethod(decompDict),
    coeffsDict_(decompDict.optionalSubDict(derivedType + "Coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can only use this decomposition method for entire mesh" << nl
            << "and supply one coordinate (cellCentre) for every cell." << nl
            << "The number of coordinates " << points.size() << nl
            << "The number of cells in the mesh " << mesh.nCells() << nl
            << exit(FatalError);
    }

    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    // Decompose using default weights
    List<int> decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), pointWeights, decomp);

    return labelList(decomp);
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
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

    // Decompose using default weights
    List<int> decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), agglomWeights, decomp);


    // Rework back into decomposition for original mesh
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
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

    // Decompose using default weights
    List<int> decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), cellWeights, decomp);

    return labelList(decomp);
}

// ************************************************************************* //
