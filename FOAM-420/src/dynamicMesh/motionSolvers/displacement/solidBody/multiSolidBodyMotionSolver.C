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

#include "motionSolvers/displacement/solidBody/multiSolidBodyMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/Fields/transformField/transformField.H"
#include "meshes/polyMesh/zones/ZoneMesh/cellZoneMesh.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::multiSolidBodyMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointsMotionSolver(mesh, dict, typeName)
{
    wordList cellZoneNames;
    if (!coeffDict().found("referenceFrames"))
    {
        forAllConstIter(dictionary, coeffDict(), iter)
        {
            if (iter().isDict())
            {
                cellZoneNames.append(iter().keyword());
            }
        }
    }
    else
    {
        cellZoneNames =
            coeffDict().lookup<wordList>("cellZones");
    }

    const label nZones = cellZoneNames.size();
    zoneIDs_.setSize(nZones);
    pointIDs_.setSize(nZones);

    forAll(cellZoneNames, namei)
    {
        const word& cellZoneName = cellZoneNames[namei];
        zoneIDs_[namei] = mesh.cellZones().findZoneID(cellZoneName);

        if (zoneIDs_[namei] == -1)
        {
            FatalIOErrorInFunction
            (
                coeffDict()
            )   << "Cannot find cellZone named " << cellZoneName
                << ". Valid zones are " << mesh.cellZones().names()
                << exit(FatalIOError);
        }

        // Collect points of cell zone.
        const cellZone& cellZoneNameI = mesh.cellZones()[zoneIDs_[namei]];

        boolList movePts(mesh.nPoints(), false);

        forAll(cellZoneNameI, i)
        {
            label celli = cellZoneNameI[i];
            const cell& cFaces = mesh.cells()[celli];
            forAll(cFaces, j)
            {
                const face& f = mesh.faces()[cFaces[j]];
                forAll(f, k)
                {
                    const label pointi = f[k];
                    movePts[pointi] = true;
                }
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(mesh.nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        pointIDs_[namei].transfer(ptIDs);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::~multiSolidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::multiSolidBodyMotionSolver::curPoints() const
{
    tmp<pointField> ttransformedPts(new pointField(mesh().points()));
    pointField& transformedPts = ttransformedPts.ref();

    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];

        UIndirectList<point>(transformedPts, zonePoints) =
            transformPoints
            (
                transformation(i),
                pointField(points0_, zonePoints)
            );
    }

    return ttransformedPts;
}


void Foam::multiSolidBodyMotionSolver::solve()
{
    // Support for topology change isn't included here because it will be
    // deprecated in the future by movers and topoChangers.
    if
    (
        coeffDict().lookupOrDefault<Switch>("resetPoints0", false)
     || (!coeffDict().found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "Resetting points0" << endl;
        points0() = mesh().points();
    }
}


// ************************************************************************* //
