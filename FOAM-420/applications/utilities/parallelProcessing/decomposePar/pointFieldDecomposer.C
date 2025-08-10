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

#include "pointFieldDecomposer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const pointPatch& completeMeshPatch,
    const pointPatch& procMeshPatch,
    const labelList& directAddr
)
:
    pointPatchFieldMapperPatchRef
    (
        completeMeshPatch,
        procMeshPatch
    ),
    directAddressing_(procMeshPatch.size(), -1),
    hasUnmapped_(false)
{
    // Create the inverse-addressing of the patch point labels.
    labelList pointMap(completeMeshPatch.boundaryMesh().mesh().size(), -1);

    const labelList& completeMeshPatchPoints = completeMeshPatch.meshPoints();

    forAll(completeMeshPatchPoints, pointi)
    {
        pointMap[completeMeshPatchPoints[pointi]] = pointi;
    }

    // Use the inverse point addressing to create the addressing table for this
    // patch
    const labelList& procMeshPatchPoints = procMeshPatch.meshPoints();

    forAll(procMeshPatchPoints, pointi)
    {
        directAddressing_[pointi] =
            pointMap[directAddr[procMeshPatchPoints[pointi]]];
    }

    // Check that all the patch point addresses are set
    if (directAddressing_.size() && min(directAddressing_) < 0)
    {
        hasUnmapped_ = true;

        FatalErrorInFunction
            << "Incomplete patch point addressing"
            << abort(FatalError);
    }
}


Foam::pointFieldDecomposer::pointFieldDecomposer
(
    const pointMesh& completeMesh,
    const pointMesh& procMesh,
    const labelList& pointAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    pointAddressing_(pointAddressing),
    patchFieldDecomposers_(procMesh_.boundary().size())
{
    forAll(procMesh_.boundary(), patchi)
    {
        if (patchi < completeMesh_.boundary().size())
        {
            patchFieldDecomposers_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    completeMesh_.boundary()[patchi],
                    procMesh_.boundary()[patchi],
                    pointAddressing_
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::~pointFieldDecomposer()
{}


// ************************************************************************* //
