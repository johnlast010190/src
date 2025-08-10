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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/solidBody/multiSolidBodyGIBMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/Fields/transformField/transformField.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyGIBMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyGIBMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyGIBMotionSolver::multiSolidBodyGIBMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    triName_(coeffDict().lookup("triSurfaceName", true))
{
    wordList names;

    forAllConstIter(dictionary, coeffDict(), iter)
    {
        if (iter().isDict())
        {
            names.append(iter().keyword());
        }
    }
    const wordList regionNames =
        coeffDict().lookupOrDefault<wordList>("regionNames", names);

    triSurfaceMesh stl
    (
        IOobject
        (
            triName_,
            mesh.time().constant(),
            "triSurface",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    triPoints0_ = stl.points();

    const label nZones = regionNames.size();
    zoneIDs_.setSize(nZones);
    pointIDs_.setSize(nZones);

    forAll(regionNames, zonei)
    {
        zoneIDs_[zonei] = -1;
        forAll(stl.regions(), i)
        {
            if (stl.regions()[i] == regionNames[zonei])
            {
                zoneIDs_[zonei] = i;
            }
        }

        if (zoneIDs_[zonei] == -1)
        {
            FatalIOErrorInFunction
            (
                coeffDict()
            )   << "Cannot find region named " << regionNames[zonei]
                << ". in stl " << triName_
                << exit(FatalIOError);
        }
    }

    forAll(zoneIDs_, zi)
    {
        DynamicList<label> dynPointIDs;
        forAll(stl.surfFaces(), trii)
        {
            if (zoneIDs_[zi] == stl.surfFaces()[trii].region())
            {
                forAll(stl.surfFaces()[trii], triPti)
                {
                    dynPointIDs.append(stl.surfFaces()[trii][triPti]);
                }
            }
        }
        pointIDs_[zi] = labelList(dynPointIDs);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyGIBMotionSolver::~multiSolidBodyGIBMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::multiSolidBodyGIBMotionSolver::curPoints() const
{
    tmp<pointField> ttransformedPts(new pointField(triPoints0_));
    pointField& transformedPts = ttransformedPts.ref();

    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];
        const septernion transform = transformation(i);

        UIndirectList<point>(transformedPts, zonePoints) =
            transformPoints
            (
                transform,
                pointField(triPoints0_, zonePoints)
            );
    }

    return ttransformedPts;
}


// ************************************************************************* //
