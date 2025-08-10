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

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvFieldDecomposer::completePatchID
(
    const label procPatchi
) const
{
    const label nDirectPatches =
        completeMesh_.boundary().size()
      - completeMesh_.boundaryMesh().nIndirectPatches();

    const fvPatch& procPatch = procMesh_.boundary()[procPatchi];

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

Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressing
)
:
    labelList(mag(addressing) - 1),
    directFvPatchFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::fvFieldDecomposer::indirectPatchFieldDecomposer::
indirectPatchFieldDecomposer
(
    const labelUList& procPatchAddressing,
    const labelUList& completePatchAddressing,
    const labelList& faceAddressing
)
:
    indirectAddressing_(procPatchAddressing.size(), 0)
{
    labelList gProcAddressing(procPatchAddressing.size(), 0);
    forAll(gProcAddressing, fI)
    {
        gProcAddressing[fI] = faceAddressing[procPatchAddressing[fI]] - 1;
    }

    forAll(indirectAddressing_, fI)
    {
        forAll(completePatchAddressing, cfI)
        {
            if (gProcAddressing[fI] == completePatchAddressing[cfI])
            {
                indirectAddressing_[fI] = cfI;
            }
        }
    }
}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const surfaceLabelField::Boundary& faceAddressingBf
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    faceAddressingBf_(faceAddressingBf),
    patchFieldDecomposers_(procMesh_.boundary().size()),
    indirectPatchFieldDecomposers_(procMesh_.boundary().size())
{
    forAll(procMesh_.boundary(), procPatchi)
    {
        const fvPatch& procPatch = procMesh_.boundary()[procPatchi];

        if (isA<directPolyPatch>(procPatch.patch()))
        {
            // Determine the index of the corresponding complete patch
            const label completePatchi = completePatchID(procPatchi);

            // If there is a corresponding complete patch, then create
            // a patch mapper
            if (completePatchi >= 0)
            {
                patchFieldDecomposers_.set
                (
                    procPatchi,
                    new patchFieldDecomposer(faceAddressingBf[completePatchi])
                );
            }
        }
        else
        {
            const polyPatch& pp = procMesh_.boundary()[procPatchi].patch();
            const indirectPolyPatch& procIndPP =
                refCast<const indirectPolyPatch>(pp);

            label completeIndID =
                completeMesh_.boundaryMesh().findPatchID(pp.name());
            const indirectPolyPatch& completeIndPP =
                refCast<const indirectPolyPatch>
                (
                    completeMesh_.boundaryMesh()[completeIndID]
                );

            indirectPatchFieldDecomposers_.set
            (
                procPatchi,
                new indirectPatchFieldDecomposer
                (
                    procIndPP.fAddr(),
                    completeIndPP.fAddr(),
                    faceAddressing_
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::~fvFieldDecomposer()
{}


// ************************************************************************* //
