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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "cellCellStencil/cellCellStencil/cellCellStencil.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCellStencil, 0);
    defineRunTimeSelectionTable(cellCellStencil, mesh);

    template<>
    const char* NamedEnum
    <
        cellCellStencil::cellType,
        3
    >::names[] =
    {
        "calculated",
        "interpolated",
        "hole"
    };
}

const Foam::NamedEnum<Foam::cellCellStencil::cellType, 3>
Foam::cellCellStencil::cellTypeNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencil::cellCellStencil(const fvMesh& mesh)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::cellCellStencil> Foam::cellCellStencil::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing cellCellStencil" << endl;
    }

    word type(dict.lookup("method"));

    const auto ctor = ctorTableLookup("cellCellStencil type", meshConstructorTable_(), type);
    return autoPtr<cellCellStencil>(ctor(mesh, dict, true));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencil::~cellCellStencil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelIOList& Foam::cellCellStencil::zoneID() const
{
    if (!mesh_.foundObject<labelIOList>("zoneID"))
    {
        labelIOList* zoneIDPtr = new labelIOList
        (
            IOobject
            (
                "zoneID",
                mesh_.facesInstance(),
                polyMesh::meshSubDir,
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_.nCells()
        );
        labelIOList& zoneID = *zoneIDPtr;

        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                mesh_.time().findInstance(mesh_.dbDir(), "zoneID"),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_
        );
        forAll(volZoneID, cellI)
        {
            zoneID[cellI] = label(volZoneID[cellI]);
        }

        zoneIDPtr->store();
    }
    return mesh_.lookupObject<labelIOList>("zoneID");
}


Foam::labelList Foam::cellCellStencil::count
(
    const label size,
    const labelUList& lst
)
{
    labelList count(size, 0);
    forAll(lst, i)
    {
        count[lst[i]]++;
    }
    Pstream::listCombineGather(count, plusEqOp<label>());
    return count;
}


bool Foam::cellCellStencil::localStencil(const labelUList& slots) const
{
    forAll(slots, i)
    {
        if (slots[i] >= mesh_.nCells())
        {
            return false;
        }
    }
    return true;
}


void Foam::cellCellStencil::globalCellCells
(
    const globalIndex& gi,
    const polyMesh& mesh,
    const labelList& selectedCells,
    labelListList& cellCells,
    pointListList& cellCellCentres
)
{
    // For selected cells determine the face neighbours (in global numbering)

    const pointField& cellCentres = mesh.cellCentres();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const cellList& cells = mesh.cells();


    // 1. Determine global cell number on other side of coupled patches

    labelList globalCellIDs(mesh.nCells());
    forAll(globalCellIDs, celli)
    {
        globalCellIDs[celli] = gi.toGlobal(celli);
    }

    labelList nbrGlobalCellIDs;
    syncTools::swapBoundaryCellList
    (
        mesh,
        globalCellIDs,
        nbrGlobalCellIDs
    );
    pointField nbrCellCentres;
    syncTools::swapBoundaryCellList
    (
        mesh,
        cellCentres,
        nbrCellCentres
    );


    // 2. Collect cell and all its neighbours

    cellCells.setSize(mesh.nCells());
    cellCellCentres.setSize(cellCells.size());

    forAll(selectedCells, i)
    {
        label celli = selectedCells[i];

        const cell& cFaces = cells[celli];
        labelList& stencil = cellCells[celli];
        pointList& stencilPoints = cellCellCentres[celli];
        stencil.setSize(cFaces.size()+1);
        stencilPoints.setSize(stencil.size());
        label compacti = 0;

        // First entry is cell itself
        stencil[compacti] = globalCellIDs[celli];
        stencilPoints[compacti++] = cellCentres[celli];

        // Other entries are cell neighbours
        forAll(cFaces, i)
        {
            label facei = cFaces[i];
            label bFacei = facei-mesh.nInternalFaces();
            label own = faceOwner[facei];
            label nbrCelli;
            point nbrCc;
            if (bFacei >= 0)
            {
                nbrCelli = nbrGlobalCellIDs[bFacei];
                nbrCc = nbrCellCentres[bFacei];
            }
            else
            {
                if (own != celli)
                {
                    nbrCelli = gi.toGlobal(own);
                    nbrCc = cellCentres[own];
                }
                else
                {
                    label nei = faceNeighbour[facei];
                    nbrCelli = gi.toGlobal(nei);
                    nbrCc = cellCentres[nei];
                }
            }

            SubList<label> current(stencil, compacti);
            if (findIndex(current, nbrCelli) == -1)
            {
                stencil[compacti] = nbrCelli;
                stencilPoints[compacti++] = nbrCc;
            }
        }
        stencil.setSize(compacti);
        stencilPoints.setSize(compacti);
    }
}


// ************************************************************************* //
