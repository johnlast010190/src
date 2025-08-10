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

#include "volumeVelocityLimit/volumeVelocityLimitObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeVelocityLimitObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        volumeVelocityLimitObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * *Protected Member Functions * * * * * * * * * * * //

void
Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
addCostOfCell
(
    scalar magU,
    scalar V
)
{
    if (magU > velocityLimitValue_ - SMALL)
    {
        objectiveValue_ += (magU - velocityLimitValue_)*V;
    }
}


void
Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
addCellInUaSource
(
    vector& source,
    vector U,
    scalar V
) const
{
    if (mag(U) > velocityLimitValue_ - SMALL)
    {
        source -= U/(mag(U) + SMALL)*V;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
volumeVelocityLimitObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    velocityLimitValue_
    (
        readScalar(objectiveDict.lookup("velocityLimitValue"))
    )
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
volumeVelocityLimitObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    volumeVelocityLimitObjectiveFunctionObject
    (
        name,
        runTime,
        objectiveDict,
        false
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
~volumeVelocityLimitObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    velocityLimitValue_ = readScalar(dict.lookup("velocityLimitValue"));

    return true;
}


bool
Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::
execute()
{
    const volVectorField& U = objectiveFunctionObject::U();
    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    objectiveValue_ = 0.0;

    if (zoneNames_.size() != 0)
    {
        forAll(zoneNames_, zI)
        {
            label index = mesh_.cellZones().findZoneID(zoneNames_[zI]);

            const labelList& zoneCells = mesh_.cellZones()[index];
            forAll(zoneCells, zcI)
            {
                label cI = zoneCells[zcI];

                addCostOfCell(mag(U[cI]), V[cI]);
            }
        }
    }
    else
    {
        forAll(V, cI)
        {
            addCostOfCell(mag(U[cI]), V[cI]);
        }
    }

    reduce(objectiveValue_, sumOp<scalar>());

    Info<< type() << " " << name() << " execute:" << nl
        << "Objective value = " << objectiveValue_ << " [m3]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::volumeVelocityLimitObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
