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
    (c) 2011-2015 OpenFOAM Foundation

InClass
    renumberMethod

\*---------------------------------------------------------------------------*/

#include "renumberMethod/renumberMethod.H"
#include "decompositionMethod/decompositionMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(renumberMethod, 0);
    defineRunTimeSelectionTable(renumberMethod, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::renumberMethod> Foam::renumberMethod::New
(
    const dictionary& renumberDict
)
{
    const word methodType(renumberDict.lookup("method"));

    //Info<< "Selecting renumberMethod " << methodType << endl;

    const auto ctor = ctorTableLookup("renumberMethod type", dictionaryConstructorTable_(), methodType);
    return autoPtr<renumberMethod>(ctor(renumberDict));
}


Foam::labelList Foam::renumberMethod::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    CompactListList<label> cellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        false,                      // local only
        cellCells
    );

    // Renumber based on agglomerated points
    return renumber(cellCells(), points);
}


Foam::labelList Foam::renumberMethod::renumber
(
    const labelList& cellCells,
    const labelList& offsets,
    const pointField& cc
) const
{
    NotImplemented;
}


Foam::labelList Foam::renumberMethod::renumber
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints
) const
{
    CompactListList<label> coarseCellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),
        false,                      // local only
        coarseCellCells
    );

    // Renumber based on agglomerated points
    labelList coarseDistribution
    (
        renumber
        (
            coarseCellCells(),
            coarsePoints
        )
    );

    // Rework back into renumbering for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


// ************************************************************************* //
