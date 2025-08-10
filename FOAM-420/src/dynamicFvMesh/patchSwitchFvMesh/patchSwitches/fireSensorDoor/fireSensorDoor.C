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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "patchSwitchFvMesh/patchSwitches/fireSensorDoor/fireSensorDoor.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fireSensorDoor, 0);

    addToRunTimeSelectionTable
    (
        patchSwitch,
        fireSensorDoor,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool fireSensorDoor::enableCondition()
{
    const volScalarField& ps =
        mesh_.lookupObject<volScalarField>(psName_);

    bool found = false;
    bool triggerFire = false;
    scalar cellValue(0);

    forAll(patchNames_, pI)
    {
        label lpI = mesh_.boundary().findPatchID(patchNames_[pI]);
        if (lpI>=0)
        {
            found = true;
            forAll(ps.boundaryField()[lpI], fI)
            {
                if (ps.boundaryField()[lpI][fI] > threshold_)
                {
                    triggerFire = true;
                }
                cellValue = ps.boundaryField()[lpI][fI];
                if (triggerFire) break;
            }
        }
    }
    forAll(points_, pI)
    {
        label cI = mesh_.findCell(points_[pI]);
        if (cI!=-1)
        {
            found = true;
            if (ps[cI]>threshold_)
            {
                //Pout<< ps[cI] << tab << threshold_ << endl;
                triggerFire = true;
            }
            cellValue = ps[cI];
            if (triggerFire) break;
        }
    }

    reduce(
        std::tie(found, triggerFire, cellValue),
        ParallelOp<orOp<bool>, orOp<bool>, sumOp<scalar>>{}
    );

    if (!found)
    {
        Info<< this->type()
             << " for GIB faceZone "
             << mesh_.faceZones()[zoneID_].name()
             << " Points: " << points_ << tab
             << "out of bounds" << endl;
    }
    else
    {
        Info<< this->type()
             << " for GIB faceZone "
             << mesh_.faceZones()[zoneID_].name()
             << " current value: "
             << cellValue
             << " threshold value: "
             << threshold_
             << ": door operating: "
             << found << endl;
    }

    return triggerFire;

}

bool fireSensorDoor::disableCondition()
{

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fireSensorDoor::fireSensorDoor
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    patchSwitch(mesh, dict),
    threshold_(readScalar(dict.lookup("openThreshold"))),
    patchNames_(0),
    points_(0),
    psName_(dict.lookup("fieldName"))
{
    if (dict.found("patches"))
    {
        patchNames_ = dict.lookup<wordList>("patches");
    }
    if (dict.found("locations"))
    {
        points_ = dict.lookup<List<point>>("locations");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
