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

#include "reconstructLagrangian.H"
#include "primitives/ints/lists/labelIOList.H"
#include "passiveParticle/passiveParticleCloud.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::reconstructLagrangianPositions
(
    const polyMesh& mesh,
    const word& cloudName,
    const PtrList<fvMesh>& meshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing
)
{
    passiveParticleCloud lagrangianPositions
    (
        mesh,
        cloudName,
        IDLList<passiveParticle>()
    );

    forAll(meshes, i)
    {
        const labelList& cellMap = cellProcAddressing[i];

        // faceProcAddressing not required currently.
        // const labelList& faceMap = faceProcAddressing[i];

        Cloud<passiveParticle> lpi(meshes[i], cloudName, false);

        forAllConstIter(Cloud<passiveParticle>, lpi, iter)
        {
            const passiveParticle& ppi = iter();

            // // Inverting sign if necessary and subtracting 1 from
            // // faceProcAddressing
            // label mappedTetFace = mag(faceMap[ppi.tetFace()]) - 1;

            lagrangianPositions.append
            (
                new passiveParticle
                (
                    mesh,
                    ppi.position(),
                    cellMap[ppi.cell()],
                    false
                )
            );
        }
    }

    IOPosition<Cloud<passiveParticle>>(lagrangianPositions).write();
}


// ************************************************************************* //
