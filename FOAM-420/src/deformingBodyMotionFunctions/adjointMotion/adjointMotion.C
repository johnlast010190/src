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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "adjointMotion.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace deformingBodyMotionFunctions
{
    defineTypeNameAndDebug(adjointMotion, 0);
    addToRunTimeSelectionTable
    (
        deformingBodyMotionFunction,
        adjointMotion,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField> adjointMotion::normalVelocity(const label& pI) const
{
    tmp<scalarField > Unt
    (
        new scalarField(mesh_.boundary()[pI].size())
    );
    scalarField& Un = Unt.ref();

    volScalarField& G =
        const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("G"));

    forAll(Un, fI)
    {
        Un[fI] = G.boundaryField()[pI][fI];
        Un[fI] = min(Un[fI],0.3);
        Un[fI] = max(Un[fI],-0.3);
    }

    return Unt;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointMotion::adjointMotion
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

adjointMotion::~adjointMotion() {}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool adjointMotion::read
(
    const dictionary& DBMFCoeffs
)
{
    return true;
}

void adjointMotion::update()
{
}

Foam::tmp<vectorField> adjointMotion::boundaryVelocity
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

Foam::tmp<vectorField> adjointMotion::interfaceVelocity
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
