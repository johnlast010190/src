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
    (c) 2016-2020 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensightWrite/ensightWrite.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * //

namespace Foam
{
    // A limited selection of actions
    const Enum<topoSetSource::setAction> functionObjects::ensightWrite::actionNames
    ({
        { topoSetSource::NEW, "use" },  // Reuse NEW for "use" action name
        { topoSetSource::ADD, "add" },
        { topoSetSource::SUBSET, "subset" },
        { topoSetSource::INVERT, "invert" },
    });
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::updateSubset
(
    fvMeshSubset& subsetter
)
{
    if
    (
        selection_.empty()
    )
    {
        return false;
    }

    const fvMesh& mesh = subsetter.baseMesh();

    // Start with all cells unselected
    cellSet cellsToSelect(mesh, "subEnsightCellSet", mesh.nCells()/10+1);
    for (const entry& dEntry : selection_)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Ignoring non-dictionary entry "
                << dEntry << endl;
            continue;
        }

        const dictionary& dict = dEntry.dict();

        const auto action = actionNames.lookup("action", dict);

        // Handle manually
        if (action == topoSetSource::INVERT)
        {
            cellsToSelect.invert(mesh.nCells());
            continue;
        }

        auto source = topoSetSource::New
        (
            dict.lookup("source"),
            mesh,
            dict.optionalSubDict("sourceInfo")
        );

        switch (action)
        {
            case topoSetSource::NEW:  // "use"
            case topoSetSource::ADD:
                if (topoSetSource::NEW == action)
                {
                }
                source->applyToSet(action, cellsToSelect);
                break;

            case topoSetSource::SUBSET:
            {
                cellSet other(mesh, "writeEnsightUpdateOtherCellSet");
                source->applyToSet(topoSetSource::NEW, other);

                cellsToSelect.subset(other);
            }
            break;

            default:
                // Should already have been caught
                WarningInFunction
                    << "Ignoring unhandled action '"
                    << actionNames[action] << "'" << endl;
                break;
        }
    }

    subsetter.clear();
    subsetter.setCellSubset(static_cast<labelHashSet>(cellsToSelect));

    ensMesh_.clear();

    return true;
}


bool Foam::functionObjects::ensightWrite::update()
{
    if (meshState_ == polyMesh::UNCHANGED)
    {
        return false;
    }

    updateSubset(meshSubset_);

    meshState_ = polyMesh::UNCHANGED;

    if (!ensMesh_.valid())
    {
        ensMesh_.reset(new ensightMesh(meshSubset_.mesh(), writeOpts_));
    }
    else if
    (
        ensMesh_().needsUpdate()
    )
    {
        ensMesh_().correct();
    }

    return true;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::ensightWrite::readSelection(const dictionary& dict)
{
    // Ensure consistency
    ensMesh_.clear();

    meshSubset_.clear();
    meshState_ = polyMesh::TOPO_CHANGE;

    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh& mpm)
{
    meshState_ = polyMesh::TOPO_CHANGE;

    if (ensMesh_.valid())
    {
        ensMesh_().expire();
    }
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh& mpm)
{
    // Only move to worse states
    if (meshState_ == polyMesh::UNCHANGED)
    {
        meshState_ = polyMesh::POINTS_MOVED;
    }
    if (ensMesh_.valid())
    {
        ensMesh_().expire();
    }
}


// ************************************************************************* //
