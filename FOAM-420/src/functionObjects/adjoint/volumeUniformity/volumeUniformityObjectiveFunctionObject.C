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

#include "volumeUniformity/volumeUniformityObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeUniformityObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        volumeUniformityObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
zoneVolume() const
{
    const DimensionedField<scalar, volMesh>& V(mesh_.V());

    scalar volume = 0;

    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells
        (
            mesh_.cellZones()[index]
        );

        forAll(sCells, sI)
        {
            label i = sCells[sI];
            volume  += V[i];
        }
    }

    reduce(volume, sumOp<scalar>());

    if (volume == 0)
    {
        FatalErrorInFunction
            << "Volume uniformity zone contains no cells."
            << exit(FatalError);
    }

    return volume;
}


Foam::vector
Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
targetVelocity() const
{
    if (targetVelocity_.valid())
    {
        return targetVelocity_();
    }
    else
    {
        const volVectorField& u(U());
        const DimensionedField<scalar, volMesh>& V(mesh_.V());

        scalar UmeanDir = 0;

        forAll(zoneNames_, szI)
        {
            label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
            const labelList& sCells
            (
                mesh_.cellZones()[index]
            );

            forAll(sCells, sI)
            {
                label i = sCells[sI];
                UmeanDir  += V[i]*(u[i] & desiredFlowDirection_());
            }
        }

        reduce(UmeanDir, sumOp<scalar>());
        UmeanDir /= volume_;

        return (UmeanDir * desiredFlowDirection_());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
volumeUniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    desiredFlowDirection_(),
    targetVelocity_(),
    sourceCoeff_(objectiveDict.lookupOrDefault<scalar>("powerScale", 1)),
    volume_(zoneVolume())
{
    createFiles(useAdjointFileFormat);

    if (objectiveDict.found("targetVelocity"))
    {
        targetVelocity_.reset
        (
            new vector(objectiveDict.lookup("targetVelocity"))
        );
    }
    else
    {
        desiredFlowDirection_.reset
        (
            new vector(objectiveDict.lookup("uniformityDirection"))
        );
        desiredFlowDirection_() /= mag(desiredFlowDirection_());
    }
}


Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
volumeUniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    volumeUniformityObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
~volumeUniformityObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::volumeUniformityObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    sourceCoeff_ = dict.lookupOrDefault<scalar>("powerScale", 1);
    volume_ = zoneVolume();

    if (dict.found("targetVelocity"))
    {
        targetVelocity_.reset(new vector(dict.lookup("targetVelocity")));
    }
    else
    {
        desiredFlowDirection_.reset
        (
            new vector(dict.lookup("uniformityDirection"))
        );
        desiredFlowDirection_() /= mag(desiredFlowDirection_());
    }

    return true;
}


bool
Foam::functionObjects::volumeUniformityObjectiveFunctionObject::
execute()
{
    const volVectorField& V(U());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());

    objectiveValue_ = 0;
    vector Vtarget = targetVelocity();
    scalar meanDeviation = 0;

    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells
        (
            mesh_.cellZones()[index]
        );

        forAll(sCells, sI)
        {
            label i = sCells[sI];

            objectiveValue_ += 0.5*Vol[i]*(sourceCoeff_*magSqr(V[i] - Vtarget)
                +(1-sourceCoeff_)*sqr(magSqr(V[i] - Vtarget)));

            meanDeviation += Vol[i]*mag(V[i] - Vtarget);
        }
    }

    reduce(objectiveValue_, sumOp<scalar>());
    reduce(meanDeviation, sumOp<scalar>());

    meanDeviation /= volume_;

    Info<< type() << " " << name() << " execute:" << nl
        << "Velocity deviation = " << objectiveValue_
        << " [m/s]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::volumeUniformityObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
