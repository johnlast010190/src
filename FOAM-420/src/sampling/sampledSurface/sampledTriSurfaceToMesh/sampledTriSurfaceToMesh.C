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
    (c) 2017 OpenCFD Ltd.
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledTriSurfaceToMesh/sampledTriSurfaceToMesh.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "db/dictionary/dictionary.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledTriSurfaceToMesh, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledTriSurfaceToMesh,
        word,
        triSurfaceToMesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceToMesh::sampledTriSurfaceToMesh
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const keyType& zoneKey,
    const bool planarCut,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingTriSurface
    (
        triSurfaceMesh
        (
            IOobject
            (
                surfaceName,
                mesh.time().constant(), // instance
                "triSurface",           // local
                mesh,                   // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
             )
        )
    ),
    fvMesh_(refCast<const fvMesh>(mesh)),
    coorFramePtr_(&coordinateFrame::New(fvMesh_, name)),
    zoneKey_(zoneKey),
    triangulate_(triangulate),
    planarCut_(planarCut),
    needsUpdate_(true),
    moveSurface_(false),
    motionPoints_(surface_.triSurface::points())
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }

    // Try to read point values when simulation is re-started
    if (mesh.time().value() != 0.0 && coorFramePtr_->isIncrementalMotion())
    {
        IOobject header
        (
            name + "Points",
            mesh.time().timeName(),
            "uniform",
            mesh.parent(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );
        if (header.typeHeaderOk<IOdictionary>(true))
        {
            IOdictionary oldPointsDict(header);
            oldPointsDict.lookup("surfacePoints") >> motionPoints_;
        }

        coorFramePtr_->updateState();

        const pointField newPoints
        (
            transformPoints
            (
                coorFramePtr_->transformation(),
                surface_.triSurface::points()
            )
        );

        surface_.triSurface::movePoints(newPoints);
    }
}


Foam::sampledTriSurfaceToMesh::sampledTriSurfaceToMesh
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingTriSurface
    (
        triSurfaceMesh
        (
            IOobject
            (
                dict.lookup("surface"),
                mesh.time().constant(), // instance
                "triSurface",           // local
                mesh,                   // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    ),
    fvMesh_(refCast<const fvMesh>(mesh)),
    coorFramePtr_
    (
        dict.found("referenceFrame")
      ? coordinateFrame::lookupNew(fvMesh_, dict)
      : nullptr
    ),
    zoneKey_(keyType::null),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    planarCut_(dict.lookupOrDefault("planarCutMethod", false)),
    needsUpdate_(true),
    moveSurface_(dict.lookupOrDefault<Switch>("moveSurface", false)),
    motionPoints_(surface_.triSurface::points())
{
    dict.readIfPresent("zone", zoneKey_);

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }

    if (coorFramePtr_)
    {
        //coorFramePtr_ = coordinateFrame::lookupNew(fvMesh_, dict);
        // Try to read point values when simulation is re-started
        if (mesh.time().value() != 0.0 && coorFramePtr_->isIncrementalMotion())
        {
            IOobject header
            (
                name + "Points",
                mesh.time().timeName(),
                "uniform",
                mesh.parent(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            //Info<< "header " << header << endl;
            if (header.typeHeaderOk<IOdictionary>(true))
            {
                IOdictionary oldPointsDict(header);
                dictionary(oldPointsDict).lookup("surfacePoints") >> motionPoints_; //surface_.triSurface::points();
            }

            coorFramePtr_->updateState();

            const pointField newPoints
            (
                transformPoints
                (
                    coorFramePtr_->transformation(),
                    motionPoints_
                )
            );

            surface_.triSurface::movePoints(newPoints);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceToMesh::~sampledTriSurfaceToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledTriSurfaceToMesh::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledTriSurfaceToMesh::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledTriSurfaceToMesh::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    if (coorFramePtr_)
    {
        coorFramePtr_->updateState();

        const pointField newPoints
        (
            transformPoints
            (
                coorFramePtr_->transformation(),
                surface_.triSurface::points()
            )
        );

        surface_.triSurface::movePoints(newPoints);

        if (coorFramePtr_->isIncrementalMotion())
        {
            write();
        }
    }

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();

    if (returnReduce(selectedCells.empty(), andOp<bool>()))
    {
        reCut(mesh(), triangulate_, planarCut_);
    }
    else
    {
        reCut(mesh(), triangulate_, planarCut_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    if (moveSurface_)
    {
        needsUpdate_ = true;
    }
    else
    {
        needsUpdate_ = false;
    }

    return true;
}


bool Foam::sampledTriSurfaceToMesh::write()
{
    if
    (
        mesh().time().writeTime()
     && coorFramePtr_
     && coorFramePtr_->isIncrementalMotion()
    )
    {
        IOdictionary dict
        (
            IOobject
            (
                name() + "Points",
                mesh().time().timeName(),
                "uniform",
                mesh().parent(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        dict.set("surfacePoints", surface_.triSurface::points());
        dict.regIOobject::write();
    }

    return true;
}

Foam::tmp<Foam::scalarField> Foam::sampledTriSurfaceToMesh::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledTriSurfaceToMesh::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledTriSurfaceToMesh::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledTriSurfaceToMesh::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledTriSurfaceToMesh::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledTriSurfaceToMesh::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledTriSurfaceToMesh::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledTriSurfaceToMesh::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledTriSurfaceToMesh::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledTriSurfaceToMesh::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledTriSurfaceToMesh::print(Ostream& os) const
{
    os  << "sampledTriSurfaceToMesh: " << name() << " :"
        << "  triangulate:" << triangulate_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
