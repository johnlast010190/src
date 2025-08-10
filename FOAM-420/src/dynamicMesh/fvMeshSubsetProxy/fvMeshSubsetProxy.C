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
    (c) 2016-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubsetProxy/fvMeshSubsetProxy.H"
#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubsetProxy::fvMeshSubsetProxy(fvMesh& baseMesh)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(-1),
    type_(NONE),
    name_(),
    names_(),
    selectedCells_()
{
    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const subsetType type,
    const word& selectionName,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(selectionName.empty() ? NONE : type),
    name_(),
    names_(),
    selectedCells_()
{
    if (type_ == ZONES)
    {
        // Populate wordRes for ZONES
        names_.resize(1);
        names_.first() = selectionName;
    }
    else if (type_ != NONE)
    {
        name_ = selectionName;
    }

    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const wordRes& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(ZONES),
    name_(),
    names_(zoneNames),
    selectedCells_()
{
    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    wordRes&& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(ZONES),
    name_(),
    names_(std::move(zoneNames)),
    selectedCells_()
{
    if (useSubMesh())
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshSubsetProxy::correct(bool verbose)
{
    if (type_ == NONE)
    {
        subsetter_.clear();
        selectedCells_.clearStorage();
        return false;
    }

    const label nCells = baseMesh_.nCells();

    PackedBoolList selectedCells;

    if (type_ == SET)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet cset(baseMesh_, name_);

        selectedCells.resize(nCells);
        for (const label idx : cset)
        {
            selectedCells.set(idx);
        }
    }
    else if (type_ == ZONE)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        selectedCells.resize(nCells);
        selectedCells.set(baseMesh_.cellZones()[name_]);
    }
    else if (type_ == ZONES)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZones "
                << flatOutput(names_) << endl;
        }

        selectedCells = baseMesh_.cellZones().findMatching(names_);
    }

    const bool changed = (selectedCells_ != selectedCells);

    // Use as a cached value for next time
    selectedCells_.transfer(selectedCells);

    if (changed || selectedCells_.empty())
    {
        // 1906 implements setCellSubset(PackedBoolList, label)
        // This has not yet been backported, so we must work around it
        // subsetter_.setCellSubset(selectedCells_, exposedPatchId_);

        // Construct a labelHashSet of all indices of selectedCells_ where
        // the value is "true"
        labelHashSet temporaryHashSet;
        forAll(selectedCells_, i)
        {
            bool bool_value = selectedCells_[i];
            if (bool_value)
            {
                temporaryHashSet.insert(i);
            }
        }
        subsetter_.setCellSubset(temporaryHashSet, exposedPatchId_);
    }

    return returnReduce(changed, orOp<bool>());
}


Foam::polyMesh::readUpdateState Foam::fvMeshSubsetProxy::readUpdate()
{
    const polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if (meshState == polyMesh::POINTS_MOVED)
    {
        if (correct(true))
        {
            // The cellSet/cellZone changed on POINTS_MOVED,
            // treat like TOPO_CHANGE
            return polyMesh::TOPO_CHANGE;
        }
    }
    else if
    (
        meshState == polyMesh::TOPO_CHANGE
     || meshState == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        correct(true);
    }

    return meshState;
}


// ************************************************************************* //
