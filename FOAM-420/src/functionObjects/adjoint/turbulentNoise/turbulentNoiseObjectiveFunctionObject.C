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

#include "turbulentNoise/turbulentNoiseObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulentNoiseObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        turbulentNoiseObjectiveFunctionObject,
        dictionary
    );
}
}

// * * * * * * * * * * * * * *Protected Member Functions * * * * * * * * * * //

Foam::labelList
Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::
findZoneIDs()
{
    labelList aeroAcousticsZoneLabels(aeroAcousticsZoneNames_.size());

    forAll(aeroAcousticsZoneNames_, zoneI)
    {
        label cellZoneID
        (
            mesh_.cellZones().findZoneID
            (
                aeroAcousticsZoneNames_[zoneI]
            )
        );

        if (cellZoneID == -1 && !Pstream::parRun())
        {
            WarningInFunction
                << "cannot find aeroacoustics cellZone "
                << aeroAcousticsZoneNames_[zoneI];
        }
        aeroAcousticsZoneLabels[zoneI] = cellZoneID;
    }
    return aeroAcousticsZoneLabels;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::
turbulentNoiseObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    aeroAcousticsZoneNames_ (objectiveDict.lookup("aeroAcousticsZones")),
    aeroAcousticsZoneLabels_(findZoneIDs())
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::
turbulentNoiseObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    turbulentNoiseObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::
~turbulentNoiseObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    aeroAcousticsZoneNames_  = wordList(dict.lookup("aeroAcousticsZones"));
    aeroAcousticsZoneLabels_ = findZoneIDs();

    return true;
}


bool
Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::
execute()
{
    objectiveValue_ = 0.0;

    const volScalarField mut( muEff() - mu() );
    const volScalarField nut( mut/rho() );

    forAll(aeroAcousticsZoneLabels_, zoneI)
    {
        const label& cellZoneID   = aeroAcousticsZoneLabels_[zoneI];
        const cellZone& zoneCells = mesh_.cellZones()[cellZoneID];

        forAll(zoneCells, zcI)
        {
            const label& cellI = zoneCells[zcI];
            objectiveValue_   += sqr(nut[cellI])*(mesh_.V()[cellI]);
        }
    }
    reduce(objectiveValue_, sumOp<scalar>());

    Info<< type() << " " << name() << " execute:" << nl
        << "Turbulent noise = " << objectiveValue_ << " [m7s-2]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::turbulentNoiseObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
