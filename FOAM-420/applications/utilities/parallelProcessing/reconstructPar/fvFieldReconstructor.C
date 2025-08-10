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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvFieldReconstructor.H"
#include "fvMesh/fvPatches/constraint/processorCyclic/processorCyclicFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvFieldReconstructor::completePatchID
(
    const label proci,
    const label procPatchi
) const
{
    const label nDirectPatches =
        completeMesh_.boundary().size()
      - completeMesh_.boundaryMesh().nIndirectPatches();

    const fvPatch& procPatch = procMeshes_[proci].boundary()[procPatchi];

    if (procPatchi < nDirectPatches)
    {
        return procPatchi;
    }
    else if (isA<processorCyclicFvPatch>(procPatch))
    {
        return refCast<const processorCyclicFvPatch>(procPatch).referPatchID();
    }
    else
    {
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldReconstructor::fvFieldReconstructor
(
    const fvMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing),
    faceProcAddressingBf_(faceProcAddressingBf),
    nReconstructed_(0)
{
    forAll(procMeshes_, proci)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        if
        (
            faceProcAddressing[proci].size() != procMesh.nFaces()
         || cellProcAddressing[proci].size() != procMesh.nCells()
        )
        {
            FatalErrorInFunction
                << "Size of maps does not correspond to size of mesh"
                << " for processor " << proci << endl
                << "faceProcAddressing : " << faceProcAddressing[proci].size()
                << " nFaces : " << procMesh.nFaces() << endl
                << "cellProcAddressing : " << cellProcAddressing[proci].size()
                << " nCell : " << procMesh.nCells() << endl
                << exit(FatalError);
        }
    }
}


// ************************************************************************* //
