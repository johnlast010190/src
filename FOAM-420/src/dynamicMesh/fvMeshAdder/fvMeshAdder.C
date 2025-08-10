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

#include "fvMesh/fvMesh.H"
#include "fvMeshAdder/fvMeshAdder.H"
#include "polyMeshAdder/faceCoupleInfo.H"
#include "fvMesh/fvMesh.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(fvMeshAdder, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fvMeshAdder::calcPatchMap
(
    const label oldStart,
    const label oldSize,
    const labelList& oldToNew,
    const polyPatch& newPatch,
    const label unmappedValue
)
{
    labelList newToOld(newPatch.size(), unmappedValue);

    label newStart = newPatch.start();
    label newSize = newPatch.size();

    for (label i = 0; i < oldSize; i++)
    {
        label newFacei = oldToNew[oldStart+i];

        if (newFacei >= newStart && newFacei < newStart+newSize)
        {
            newToOld[newFacei-newStart] = i;
        }
    }
    return newToOld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::fvMeshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary,
    const bool fullyMapped
)
{
    mesh0.clearOut();

    // Resulting merged mesh (polyMesh only!)
    autoPtr<mapAddedPolyMesh> mapPtr
    (
        polyMeshAdder::add
        (
            mesh0,
            mesh1,
            coupleInfo,
            validBoundary
        )
    );

    // Adjust the fvMesh part.
    const polyBoundaryMesh& patches = mesh0.boundaryMesh();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh0.boundary());
    fvPatches.setSize(patches.size());
    forAll(patches, patchi)
    {
        fvPatches.set(patchi, fvPatch::New(patches[patchi], fvPatches));
    }

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvMeshAdder::MapVolFields<scalar>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<vector>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<sphericalTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<symmTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<tensor>(mapPtr, mesh0, mesh1, fullyMapped);

    fvMeshAdder::MapVolFields<vector1>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<vector4>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<tensor4>(mapPtr, mesh0, mesh1, fullyMapped);


    fvMeshAdder::MapSurfaceFields<scalar>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<vector>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<sphericalTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<symmTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<tensor>(mapPtr, mesh0, mesh1, fullyMapped);


    fvMeshAdder::MapDimFields<scalar>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<vector>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<symmTensor>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<tensor>(mapPtr, mesh0, mesh1);

    fvMeshAdder::MapDimFields<vector1>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<vector4>(mapPtr, mesh0, mesh1);
    fvMeshAdder::MapDimFields<tensor4>(mapPtr, mesh0, mesh1);

    return mapPtr;
}


// ************************************************************************* //
