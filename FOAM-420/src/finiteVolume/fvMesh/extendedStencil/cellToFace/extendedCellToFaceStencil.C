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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencil.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "containers/Lists/SortableList/SortableList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(extendedCellToFaceStencil, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::extendedCellToFaceStencil::writeStencilStats
(
    Ostream& os,
    const labelListList& stencil,
    const mapDistribute& map
)
{
    label sumSize = 0;
    label nSum = 0;
    label minSize = labelMax;
    label maxSize = labelMin;

    forAll(stencil, i)
    {
        const labelList& sCells = stencil[i];

        if (sCells.size() > 0)
        {
            sumSize += sCells.size();
            nSum++;
            minSize = min(minSize, sCells.size());
            maxSize = max(maxSize, sCells.size());
        }
    }
    reduce(
        std::tie(sumSize, nSum, minSize, maxSize),
        ParallelOp<sumOp<label>, sumOp<label>, minOp<label>, maxOp<label>>{}
    );

    os  << "Stencil size :" << nl
        << "    average : " << scalar(sumSize)/nSum << nl
        << "    min     : " << minSize << nl
        << "    max     : " << maxSize << nl
        << endl;

    // Sum all sent data
    label nSent = 0;
    label nLocal = 0;
    forAll(map.subMap(), proci)
    {
        if (proci != Pstream::myProcNo())
        {
            nSent += map.subMap()[proci].size();
        }
        else
        {
            nLocal += map.subMap()[proci].size();
        }
    }

    os  << "Local data size : " << returnReduceToMaster(nLocal, sumOp<label>()) << nl
        << "Sent data size  : " << returnReduceToMaster(nSent, sumOp<label>()) << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCellToFaceStencil::extendedCellToFaceStencil(const polyMesh& mesh)
:
    mesh_(mesh)
{
    // Check for transformation - not supported.
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& cpp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            if (cpp.transform().transformsPosition())
            {
                FatalErrorInFunction
                    << "Coupled patches with transformations not supported."
                    << endl
                    << "Problematic patch " << cpp.name() << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
