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
    (c) 2011-2014 OpenFOAM Foundation
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.
    (c) 2016-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh/multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh.H"

#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/Fields/transformField/transformField.H"
#include "meshes/polyMesh/zones/ZoneMesh/cellZoneMesh.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#define  DEBUG Info<<"found  in file "<<__FILE__<<", at line "<<__LINE__<<endl;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh::
multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh
(
    const IOobject& io
)
:
    dynamicMultiCriterionRefineFvMesh(io),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            io.time().findInstance(this->meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    )
{
    IOdictionary dynDict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& dynamicMeshCoeffs = dynDict.subDict(typeName + "Coeffs");

    if (undisplacedPoints_.size() != nPoints())
    {
        FatalErrorInFunction
            << dynamicMeshCoeffs
            << "Read " << undisplacedPoints_.size()
            << " undisplaced points from " << undisplacedPoints_.objectPath()
            << " but the current mesh has " << nPoints()
            << exit(FatalIOError);
    }


    zoneIDs_.setSize(dynamicMeshCoeffs.size());
    motionPtr_.setSize(dynamicMeshCoeffs.size());
    pointIDs_.setSize(dynamicMeshCoeffs.size());
    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            word zoneName(subDict.lookup("cellZone"));

            zoneIDs_[zoneI] = cellZones().findZoneID(zoneName);

            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorInFunction
                (
                    dynamicMeshCoeffs
                )   << "Cannot find cellZone named " << zoneName
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            IOobject io(dynDict);
            io.readOpt() = IOobject::NO_READ;

            motionPtr_.set(zoneI, motionSolver::New(*this, subDict));

            // Collect points of cell zone.
            const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

            boolList movePts(nPoints(), false);

            forAll(cz, i)
            {
                label cellI = cz[i];
                const cell& c = cells()[cellI];
                forAll(c, j)
                {
                    const face& f = faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointI = f[k];
                        movePts[pointI] = true;
                    }
                }
            }

            syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying motionSolver " << motionPtr_[zoneI].type()
                << " to "
                << returnReduceToMaster(pointIDs_[zoneI].size(), sumOp<label>())
                << " points of cellZone " << zoneName << endl;

            zoneI++;
        }
    }
    zoneIDs_.setSize(zoneI);
    motionPtr_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh::
~multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh::update()
{
    // refine and balance first
    dynamicMultiCriterionRefineFvMesh::update();

    static bool hasWarned = false;

    // nPoints may change so an update needed
    if (this->changing())
    {
        undisplacedPoints_ = this->points();
    }

    pointField transformedPts(undisplacedPoints_);

    // update the point IDs if the mesh has changed
    if (this->changing())
    {
        this->setV0() = this->V();
        updatePointIDs();
    }

    forAll(motionPtr_, zoneI)
    {
        tmp<pointField> tnewPoints(motionPtr_[zoneI].newPoints());
        const pointField& newPoints = tnewPoints();

        const labelList& zonePoints = pointIDs_[zoneI];
        forAll(zonePoints, i)
        {
            label pointI = zonePoints[i];
            transformedPts[pointI] = newPoints[pointI];
        }
    }

    fvMesh::movePoints(transformedPts);

    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }
    moving(true);
    topoChanging(true);

    return true;
}

void Foam::multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh::updatePointIDs()
{

    forAll(zoneIDs_,zoneI)
    {
        // Collect points of cell zone.
        const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

        boolList movePts(nPoints(), false);

        forAll(cz, i)
        {
            label cellI = cz[i];
            const cell& c = cells()[cellI];
            forAll(c, j)
            {
                const face& f = faces()[c[j]];
                forAll(f, k)
                {
                    label pointI = f[k];
                    movePts[pointI] = true;
                }
            }
        }

        syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        pointIDs_[zoneI].transfer(ptIDs);
    }
}

// ************************************************************************* //
