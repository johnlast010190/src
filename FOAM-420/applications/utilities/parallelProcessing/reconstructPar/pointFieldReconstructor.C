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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFieldReconstructor::pointFieldReconstructor
(
    const pointMesh& completeMesh,
    const PtrList<pointMesh>& procMeshes,
    const labelListList& pointProcAddressing
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    pointProcAddressing_(pointProcAddressing),
    patchPointAddressing_(procMeshes.size()),
    nReconstructed_(0)
{
    // Inverse-addressing of the patch point labels.
    labelList pointMap(completeMesh_.size(), -1);

    // Create the pointPatch addressing
    forAll(procMeshes_, proci)
    {
        const pointMesh& procMesh = procMeshes_[proci];

        patchPointAddressing_[proci].setSize(procMesh.boundary().size());

        forAll(procMesh.boundary(), patchi)
        {
            if (patchi < completeMesh_.boundary().size())
            {
                labelList& procPatchAddr =
                    patchPointAddressing_[proci][patchi];
                procPatchAddr.setSize(procMesh.boundary()[patchi].size(), -1);

                const labelList& patchPointLabels =
                    completeMesh_.boundary()[patchi].meshPoints();

                // Create the inverse-addressing of the patch point labels.
                forAll(patchPointLabels, pointi)
                {
                    pointMap[patchPointLabels[pointi]] = pointi;
                }

                const labelList& procPatchPoints =
                    procMesh.boundary()[patchi].meshPoints();

                forAll(procPatchPoints, pointi)
                {
                    procPatchAddr[pointi] =
                        pointMap
                        [
                            pointProcAddressing_[proci]
                                [procPatchPoints[pointi]]
                        ];
                }

                if (procPatchAddr.size() && min(procPatchAddr) < 0)
                {
                    FatalErrorInFunction
                        << "Incomplete patch point addressing"
                        << abort(FatalError);
                }
            }
        }
    }
}


// ************************************************************************* //
