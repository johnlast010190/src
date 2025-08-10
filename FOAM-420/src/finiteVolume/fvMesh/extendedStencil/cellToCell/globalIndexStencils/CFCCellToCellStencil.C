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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/extendedStencil/cellToCell/globalIndexStencils/CFCCellToCellStencil.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CFCCellToCellStencil::calcFaceBoundaryData
(
    labelList& neiGlobal
) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();

    neiGlobal.setSize(nBnd);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label facei = pp.start();

        if (pp.coupled())
        {
            // For coupled faces get the cell on the other side
            forAll(pp, i)
            {
                label bFacei = facei-mesh().nInternalFaces();
                neiGlobal[bFacei] = globalNumbering().toGlobal(own[facei]);
                facei++;
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label bFacei = facei-mesh().nInternalFaces();
                neiGlobal[bFacei] = -1;
                facei++;
            }
        }
        else
        {
            // For noncoupled faces get the boundary face.
            forAll(pp, i)
            {
                label bFacei = facei-mesh().nInternalFaces();
                neiGlobal[bFacei] =
                    globalNumbering().toGlobal(mesh().nCells()+bFacei);
                facei++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh(), neiGlobal);
}


void Foam::CFCCellToCellStencil::calcCellStencil
(
    labelListList& globalCellCells
) const
{
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();


    // Calculate coupled neighbour (in global numbering)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList neiGlobal(nBnd);
    calcFaceBoundaryData(neiGlobal);


    // Determine cellCells in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalCellCells.setSize(mesh().nCells());
    forAll(globalCellCells, celli)
    {
        const cell& cFaces = mesh().cells()[celli];

        labelList& cCells = globalCellCells[celli];

        cCells.setSize(cFaces.size()+1);

        label nNbr = 0;

        // Myself
        cCells[nNbr++] = globalNumbering().toGlobal(celli);

        // Collect neighbouring cells/faces
        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh().isInternalFace(facei))
            {
                label nbrCelli = own[facei];
                if (nbrCelli == celli)
                {
                    nbrCelli = nei[facei];
                }
                cCells[nNbr++] = globalNumbering().toGlobal(nbrCelli);
            }
            else
            {
                label nbrCelli = neiGlobal[facei-mesh().nInternalFaces()];
                if (nbrCelli != -1)
                {
                    cCells[nNbr++] = nbrCelli;
                }
            }
        }
        cCells.setSize(nNbr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CFCCellToCellStencil::CFCCellToCellStencil(const polyMesh& mesh)
:
    cellToCellStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    calcCellStencil(*this);
}


// ************************************************************************* //
