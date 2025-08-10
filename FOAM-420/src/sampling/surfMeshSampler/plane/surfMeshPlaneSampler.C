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

#include "surfMeshSampler/plane/surfMeshPlaneSampler.H"
#include "db/dictionary/dictionary.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fields/volFields/volFields.H"
#include "coordinate/systems/cartesianCS.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshPlaneSampler, 0);
    addNamedToRunTimeSelectionTable
    (
        surfMeshSampler,
        surfMeshPlaneSampler,
        word,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMeshPlaneSampler::transferContent()
{
    SurfaceSource& src = static_cast<SurfaceSource&>(*this);
    surfMesh& dst = getOrCreateSurfMesh();

    dst.transfer(src);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshPlaneSampler::surfMeshPlaneSampler
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const keyType& zoneKey,
    const bool triangulate
)
:
    surfMeshSampler(name, mesh),
    SurfaceSource(planeDesc),
    zoneKey_(zoneKey),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


Foam::surfMeshPlaneSampler::surfMeshPlaneSampler
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSampler(name, mesh, dict),
    SurfaceSource(plane(dict)),
    zoneKey_(keyType::null),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found(coordinateSystem::typeName_()))
    {
        coordSystem::cartesian cs
        (
            coordinateSystem::New(mesh, dict, coordinateSystem::typeName_())
        );
        plane& pln = planeDesc();

        const point  orig = cs.globalPosition(pln.refPoint());
        const vector norm = cs.globalVector(pln.normal());

        if (debug)
        {
            Info<< "plane " << name << " :"
                << " origin:" << refPoint()
                << " normal:" << normal()
                << " defined within a local coordinateSystem" << endl;
        }

        // Reassign the plane
        pln = plane(orig, norm);
    }


    dict.readIfPresent("zone", zoneKey_);

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMeshPlaneSampler::~surfMeshPlaneSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshPlaneSampler::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::surfMeshPlaneSampler::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::surfMeshPlaneSampler::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();
    if (selectedCells.empty())
    {
        reCut(mesh(), triangulate_);
    }
    else
    {
        Foam::sort(selectedCells);
        reCut(mesh(), triangulate_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    transferContent();

    needsUpdate_ = false;
    return true;
}


bool Foam::surfMeshPlaneSampler::sample
(
    const word& fieldName
) const
{
    return
    (
        sampleType<scalar>(fieldName)
     || sampleType<vector>(fieldName)
     || sampleType<sphericalTensor>(fieldName)
     || sampleType<symmTensor>(fieldName)
     || sampleType<tensor>(fieldName)
    );
}


void Foam::surfMeshPlaneSampler::print(Ostream& os) const
{
    os  << "surfMeshPlaneSampler: " << name() << " :"
        << "  base:" << cuttingPlane::refPoint()
        << "  normal:" << cuttingPlane::normal()
        << "  triangulate:" << triangulate_
        << "  faces:"  << SurfaceSource::surfFaces().size()
        << "  points:" << SurfaceSource::points().size();
}


// ************************************************************************* //
