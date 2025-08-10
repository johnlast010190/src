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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::patchDistFuncs

Description
    Collection of functions used in patch distance calculation.

SourceFiles
    patchDistFuncs.C
    patchDistFuncsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistFuncs_H
#define patchDistFuncs_H

#include "fields/Fields/Field/Field.H"
#include "fields/FieldFields/FieldField/FieldField.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyBoundaryMesh;

namespace patchDistFuncs
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
struct NoMap
{
    inline bool insert(const label, const Type&)
    {
        return false;
    }

    inline bool found(const label)
    {
        return false;
    }

    static NoMap<Type> null;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Calculate smallest true distance (and face index) from point to faces
//  wallFaces
scalar smallestDist
(
    const point& p,
    const polyPatch& patch,
    const labelUList& wallFaces,
    label& patchFacei
);

// Get point neighbours of patchFacei (including patchFacei). Note: Uses a
// linear search to determine uniqueness. For polygonal faces this might be
// quite inefficient.
void getPointNeighbours
(
    const polyPatch&,
    const label patchFacei,
    DynamicList<label>& neighbours
);

//- Size of largest patch (out of supplied subset of patches)
label maxPatchSize(const polyMesh& mesh, const labelHashSet& patchIDs);

//- Correct all cells connected to boundary (via face). Sets values in
//  wallDistCorrected. Sets nearest wall face in nearestPatchAndFace.
template<template<class> class MapType=NoMap>
void correctBoundaryFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace = NoMap<labelPair>::null
);

//- Correct all face-cells connected to boundary. Sets values in
//  wallDistCorrected. Sets nearest wall face in nearestPatchAndFace.
template<template<class> class PatchField, template<class> class MapType=NoMap>
void correctBoundaryFaceFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    FieldField<PatchField, scalar>& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace = NoMap<labelPair>::null
);

//- Correct all cells connected to wall (via point). Sets values in
//  wallDistCorrected. Uses/sets nearest wall face in nearestPatchAndFace.
template<template<class> class MapType=NoMap>
void correctBoundaryPointCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace = NoMap<labelPair>::null
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistFuncs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchDist/patchDistFuncs/patchDistFuncsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
