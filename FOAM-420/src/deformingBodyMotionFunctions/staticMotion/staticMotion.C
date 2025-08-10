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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "staticMotion.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace deformingBodyMotionFunctions
{
    defineTypeNameAndDebug(staticMotion, 0);
    addToRunTimeSelectionTable
    (
        deformingBodyMotionFunction,
        staticMotion,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


tmp<scalarField> staticMotion::normalVelocity(const label& pI) const
{
    tmp<scalarField > Unt
    (
        new scalarField(mesh_.boundary()[pI].size(), 0)
    );

    return Unt;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


staticMotion::staticMotion
(
    const fvMesh& mesh,
    const dictionary& DBMFCoeffs,
    const Time& runTime
)
:
    deformingBodyMotionFunction(mesh, DBMFCoeffs, runTime)
{
    read(DBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

staticMotion::~staticMotion() {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool staticMotion::read
(
    const dictionary& DBMFCoeffs
)
{

    return true;
}


void staticMotion::update()
{
}


Foam::tmp<vectorField> staticMotion::boundaryVelocity
(
    const label& pMaster
) const
{
    tmp<vectorField > tVelocityBoundary
    (
        new vectorField
            (
                mesh_.boundary()[pMaster].size(),
                vector::zero
            )
    );

    return tVelocityBoundary;
}


Foam::tmp<vectorField> staticMotion::interfaceVelocity
(
    const label& pMaster
) const
{
    tmp<scalarField> Unt = normalVelocity(pMaster);

    tmp<vectorField > tSpeedInter
    (
        new vectorField
            (
                Unt->size(),
                vector::zero
            )
    );

    vectorField& speedInter  = tSpeedInter.ref();

    speedInter = Unt()*mesh_.boundary()[pMaster].nf()();

    return tSpeedInter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace deformingBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
