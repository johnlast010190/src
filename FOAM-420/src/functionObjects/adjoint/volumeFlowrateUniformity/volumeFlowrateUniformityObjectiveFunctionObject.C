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

#include "volumeFlowrateUniformity/volumeFlowrateUniformityObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeFlowrateUniformityObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        volumeFlowrateUniformityObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
zoneVolume(word zoneName) const
{
    const DimensionedField<scalar, volMesh>& V(mesh_.V());

    scalar volume = 0;

    label index = mesh_.cellZones().findZoneID(zoneName);
    const labelList& sCells
    (
        mesh_.cellZones()[index]
    );

    forAll(sCells, sI)
    {
        label i = sCells[sI];
        volume  += V[i];
    }
    reduce(volume, sumOp<scalar>());

    if (volume == 0)
    {
        FatalErrorInFunction
            << "Volume massflow uniformity zone "
            << zoneName
            << " contains no cells."
            << exit(FatalError);
    }

    return volume;
}


Foam::scalar
Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
zoneThickness
(
    word zoneName,
    vector flowDir
) const
{
    scalar maxP = -GREAT;
    scalar minP = GREAT;

    label index = mesh_.cellZones().findZoneID(zoneName);
    const labelList& sCells
    (
        mesh_.cellZones()[index]
    );

    const labelListList& cellPoints = mesh_.cellPoints();
    const pointField& points = mesh_.points();

    forAll(sCells, sI)
    {
        label i = sCells[sI];

        forAll(cellPoints[i], cpI)
        {
            const point& cpoint = points[cellPoints[i][cpI]];

            scalar pDotDir(cpoint & flowDir);

            maxP = max(maxP, pDotDir);
            minP = min(minP, pDotDir);
        }
    }

    reduce(maxP, maxOp<scalar>());
    reduce(minP, minOp<scalar>());

    return (maxP - minP);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
volumeFlowrateUniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    flowDirection_(0),
    ts_(0),
    areas_(0),
    zoneTargetFraction_(0)
{
    createFiles(useAdjointFileFormat);
    read(objectiveDict);
}


Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
volumeFlowrateUniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    volumeFlowrateUniformityObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
~volumeFlowrateUniformityObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    objectivePatch_ = false;

    flowDirection_ = dict.lookupOrDefault<vectorList>("flowDirections", vectorList(0));

    if (zoneNames_.size() != flowDirection_.size())
    {
        FatalErrorInFunction
            << "Number of zones and flow directions do not match."
            << exit(FatalError);
    }

    forAll(zoneNames_, zI)
    {
        flowDirection_[zI] /= mag(flowDirection_[zI]);
    }

    // calculate derived properties
    ts_.setSize(zoneNames_.size());
    areas_.setSize(zoneNames_.size());
    zoneTargetFraction_.setSize(zoneNames_.size());
    scalar totalArea = 0;

    forAll(zoneNames_, zI)
    {
        scalar volume = zoneVolume(zoneNames_[zI]);
        ts_[zI] = zoneThickness(zoneNames_[zI], flowDirection_[zI]);
        areas_[zI] = volume / ts_[zI];
        totalArea += areas_[zI];
    }

    zoneTargetFraction_ = areas_ / totalArea;

    return true;
}


bool
Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::
execute()
{
    scalar totalInletFlowRate = inletFlowRate();

    const volVectorField& V(U());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());

    objectiveValue_ = 0.0;

    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells
        (
            mesh_.cellZones()[index]
        );

        scalar zoneFlux = 0;

        forAll(sCells, sI)
        {
            label i = sCells[sI];

            zoneFlux+= Vol[i] * (V[i] & flowDirection_[szI]);
        }

        reduce(zoneFlux, sumOp<scalar>());
        zoneFlux /= ts_[szI];

        scalar dMk = zoneFlux
            - zoneTargetFraction_[szI]*totalInletFlowRate;

        objectiveValue_ += sqr(dMk);
    }
    objectiveValue_ = sqrt(objectiveValue_);

    Info<< type() << " " << name() << " execute:" << nl
        << "Volume flow rate uniformity RMS = " << objectiveValue_
        << " [m/s]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::volumeFlowrateUniformityObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
