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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "sets/topoSets/faceSet.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(faceSet, 0);

addToRunTimeSelectionTable(topoSet, faceSet, word);
addToRunTimeSelectionTable(topoSet, faceSet, size);
addToRunTimeSelectionTable(topoSet, faceSet, set);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSet::faceSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    check(mesh.nFaces());
}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const UList<label>& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceSet::~faceSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceSet::sync(const polyMesh& mesh)
{
    boolList set(mesh.nFaces(), false);

    forAllConstIter(faceSet, *this, iter)
    {
        set[iter.key()] = true;
    }
    syncTools::syncFaceList(mesh, set, orEqOp<bool>());

    label nAdded = 0;

    forAll(set, facei)
    {
        if (set[facei])
        {
            if (insert(facei))
            {
                nAdded++;
            }
        }
        else if (found(facei))
        {
            FatalErrorInFunction
                << "Problem : syncing removed faces from set."
                << abort(FatalError);
        }
    }

    reduce(nAdded, sumOp<label>());
    if (debug && nAdded > 0)
    {
        Info<< "Added an additional " << nAdded
            << " faces on coupled patches. "
            << "(processorPolyPatch, cyclicPolyPatch)" << endl;
    }
}


Foam::label Foam::faceSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


void Foam::faceSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseFaceMap());
}


void Foam::faceSet::distribute(const mapDistributePolyMesh& map)
{
    boolList inSet(map.nOldFaces());
    forAllConstIter(faceSet, *this, iter)
    {
        inSet[iter.key()] = true;
    }
    map.distributeFaceData(inSet);

    // Count
    label n = 0;
    forAll(inSet, facei)
    {
        if (inSet[facei])
        {
            n++;
        }
    }

    clear();
    resize(n);
    forAll(inSet, facei)
    {
        if (inSet[facei])
        {
            insert(facei);
        }
    }
}


void Foam::faceSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.faceCentres(), maxLen);
}


// ************************************************************************* //
