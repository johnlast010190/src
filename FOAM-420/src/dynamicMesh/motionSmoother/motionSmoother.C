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
    (c) 2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "motionSmoother/motionSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSmoother, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const dictionary& paramDict
)
:
    motionSmootherData(pMesh),
    motionSmootherAlgo
    (
        mesh,
        pMesh,
        pp,
        motionSmootherData::displacement_,
        motionSmootherData::scale_,
        motionSmootherData::oldPoints_,
        adaptPatchIDs,
        paramDict
    )
{}


Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const pointVectorField& displacement,
    const dictionary& paramDict
)
:
    motionSmootherData(displacement),
    motionSmootherAlgo
    (
        mesh,
        const_cast<pointMesh&>(displacement.mesh()),
        pp,
        motionSmootherData::displacement_,
        motionSmootherData::scale_,
        motionSmootherData::oldPoints_,
        adaptPatchIDs,
        paramDict
    )
{}


// ************************************************************************* //
