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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/wallDist/patchDistMethods/meshWave/meshWavePatchDistMethod.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "patchDist/patchDistWave/patchDistWave.H"
#include "fvMesh/wallDist/fvPatchDistWave/fvPatchDistWave.H"
#include "fvMesh/wallDist/fvWallPointData/fvWallPointData.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(meshWave, 0);
    addToRunTimeSelectionTable(patchDistMethod, meshWave, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::meshWave::meshWave
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(dict.lookupOrDefault<Switch>("correctWalls", true))
{}


Foam::patchDistMethods::meshWave::meshWave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(correctWalls)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::meshWave::correct(volScalarField& y)
{
    y = dimensionedScalar("yWall", dimLength, GREAT);

    bool hasIndirectPatches = false;
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<indirectPolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            hasIndirectPatches = true;
            break;
        }
    }

    label nUnset;
    if (!hasIndirectPatches)
    {
        nUnset =
            fvPatchDistWave::wave<fvWallPoint>
            (
                mesh_,
                patchIDs_,
                y,
                correctWalls_
            );
    }
    else
    {
        nUnset =
            patchDistWave::wave<wallPoint>
            (
                mesh_,
                patchIDs_,
                y.primitiveFieldRef(),
                correctWalls_
            );
    }

    // Update coupled and transform BC's
    y.correctBoundaryConditions();

    return nUnset > 0;
}


bool Foam::patchDistMethods::meshWave::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    y = dimensionedScalar("yWall", dimLength, GREAT);

    const label nUnset =
        fvPatchDistWave::wave<fvWallPointData<vector>>
        (
            mesh_,
            patchIDs_,
            n.boundaryField(),
            y,
            n,
            correctWalls_
        );

    // Update coupled and transform BC's
    y.correctBoundaryConditions();
    n.correctBoundaryConditions();

    return nUnset > 0;
}


// ************************************************************************* //
