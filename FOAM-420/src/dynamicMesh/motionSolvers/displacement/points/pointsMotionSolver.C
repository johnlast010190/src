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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/points/pointsMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointsMotionSolver, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::pointsMotionSolver::points0IO(const polyMesh& mesh)
{
    const word instance =
        mesh.time().findInstance
        (
            mesh.meshDir(),
            "points0",
            IOobject::READ_IF_PRESENT
        );

    if (instance != mesh.time().constant())
    {
        // points0 written to a time folder

        return
            IOobject
            (
                "points0",
                instance,
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
    }
    else
    {
        // Check that points0 are actually in constant directory

        IOobject io
        (
            "points0",
            instance,
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (io.typeHeaderOk<pointIOField>())
        {
            return io;
        }
        else
        {
            // Copy of mesh points

            return
                IOobject
                (
                    "points",
                    mesh.time().findInstance(mesh.meshDir(), "points"),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointsMotionSolver::pointsMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(pointIOField(points0IO(mesh)))
{
    if
    (
        dict.lookupOrDefault<Switch>("resetPoints0", false)
     || (!dict.found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "resetting points0" << endl;
        points0_ = mesh.points();
    }

    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file "
            <<  typeFilePath<pointIOField>
                (
                    IOobject
                    (
                        "points",
                        mesh.time().constant(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            << exit(FatalError);
    }

    // If points0 was obtained from latest mesh points, write it out to
    // ensure stable restart behaviour.
    if (points0_.name() == "points")
    {
        fileName pointsFileName(points0_.objectPath());
        points0_.rename("points0");
        // Don't write into constant to avoid confusion
        points0_.instance() = mesh.time().timeName();
        points0_.write();
    }
}


Foam::pointsMotionSolver::pointsMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const pointIOField& points0,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(points0)
{
    if
    (
        dict.lookupOrDefault<Switch>("resetPoints0", false)
     || (!dict.found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "resetting points0" << endl;
        points0_ = mesh.points();
    }

    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file " << points0.filePath()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointsMotionSolver::~pointsMotionSolver()
{}


// ************************************************************************* //
