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
    (c) 2015 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/surfaceFields/surfaceFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class Type>
void Foam::pressurePIDControlInletVelocityFvPatchVectorField::faceZoneAverage
(
    const word& name,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    scalar& area,
    Type& average
) const
{
    const fvMesh& mesh(patch().boundaryMesh().mesh());

    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(mesh));

    const faceZone& zone = mesh.faceZones()[name];

    area = 0;
    average = Type(Zero);

    forAll(zone, faceI)
    {
        const label f(zone[faceI]);

        if (mesh.isInternalFace(f))
        {
            const scalar da(mesh.magSf()[f]);

            area += da;
            average += da*field[f];
        }
        else if (isMasterFace[f])
        {
            const label bf(f-mesh.nInternalFaces());
            const label patchID = mesh.boundaryMesh().patchID()[bf];
            const label lf(mesh.boundaryMesh()[patchID].whichFace(f));
            const scalar da(mesh.magSf().boundaryField()[patchID][lf]);

            area += da;
            average += da*field.boundaryField()[patchID][lf];
        }
    }

    reduce(area, sumOp<scalar>());
    reduce(average, sumOp<Type>());

    average /= area;
}


// ************************************************************************* //
