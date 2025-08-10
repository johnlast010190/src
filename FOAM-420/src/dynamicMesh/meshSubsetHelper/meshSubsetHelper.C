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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "meshSubsetHelper/meshSubsetHelper.H"

#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshSubsetHelper::meshSubsetHelper
(
    fvMesh& baseMesh
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    type_(NONE),
    name_()
{
    correct();
}


Foam::meshSubsetHelper::meshSubsetHelper
(
    fvMesh& baseMesh,
    const subsetType type,
    const word& name
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    type_(name.empty() ? NONE : type),
    name_(name)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshSubsetHelper::correct(bool verbose)
{
    if (type_ == SET)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet subset(baseMesh_, name_);
        subsetter_.setLargeCellSubset(subset);
    }
    else if (type_ == ZONE)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        labelHashSet subset(baseMesh_.cellZones()[name_]);
        subsetter_.setLargeCellSubset(subset);
    }
}


Foam::polyMesh::readUpdateState Foam::meshSubsetHelper::readUpdate()
{
    const polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if
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
