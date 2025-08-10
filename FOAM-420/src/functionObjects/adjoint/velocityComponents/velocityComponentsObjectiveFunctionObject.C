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

#include "velocityComponents/velocityComponentsObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(velocityComponentsObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        velocityComponentsObjectiveFunctionObject,
        dictionary
    );
}
}

// * * * * * * * * * * * * *Protected Member Functions * * * * * * * * * * * //
void Foam::functionObjects::velocityComponentsObjectiveFunctionObject::writeFileHeader
(
    Ostream& os
)
{
    writeCommented(os, "Time");
    writeDelimited(os, "axialVelocityDeviation");
    writeDelimited(os, "radialVelocityDeviation");
    writeDelimited(os, "peripheralVelocityDeviation");
    writeDelimited(os, "angleDeviation");
    writeDelimited(os, "objectiveValue");
    writeDelimited(os, "normalizedWeighting");
    os << endl;
}

bool Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
isAnActivePatch(const label& pI)
{
    const fvPatch& p = mesh_.boundary()[pI];
    return
    (
        (
            p.type() == "outlet"
         || p.patch().physicalType() == "outlet"
         || p.type() == "cyclicAMI"
        )
     && (
            objectivePatch_[pI]
        )
    );
}
void Foam::functionObjects::velocityComponentsObjectiveFunctionObject::setVariables
(
    const dictionary& objectiveDict
)
{
    if (objectiveDict.found("axialVelocity"))
    {
        dictionary dict = objectiveDict.subDict("axialVelocity");
        axialVel_.active_ = true;
        axialVel_.target_ =
            dict.lookup<scalar>("target");
        axialVel_.weights_ =
            dict.lookup<scalar>("weight");
    }

    if (objectiveDict.found("peripheralVelocity"))
    {
        dictionary dict = objectiveDict.subDict("peripheralVelocity");
        peripheralVel_.active_ = true;
        peripheralVel_.target_ =
            dict.lookup<scalar>("target");
        peripheralVel_.weights_ =
            dict.lookup<scalar>("weight");
    }

    if (objectiveDict.found("radialVelocity"))
    {
        dictionary dict = objectiveDict.subDict("radialVelocity");
        radialVel_.active_ = true;
        radialVel_.target_ =
            dict.lookup<scalar>("target");
        radialVel_.weights_ =
            dict.lookup<scalar>("weight");
    }

    if (objectiveDict.found("peripheralOverAxialVelocity"))
    {
        dictionary dict = objectiveDict.subDict("peripheralOverAxialVelocity");
        peripheralOverAxialVel_.active_ = true;
        peripheralOverAxialVel_.target_ =
            dict.lookup<scalar>("target");
        peripheralOverAxialVel_.weights_ =
            dict.lookup<scalar>("weight");
    }
}

void Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
calculateObjectiveVolume()
{
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());
    objectiveVolume_ = 0;
    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells(mesh_.cellZones()[index]);
        forAll(sCells, sI)
        {
            label i = sCells[sI];
            objectiveVolume_ += Vol[i];
        }
    }
    reduce(objectiveVolume_, sumOp<scalar>());
}

void Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
calculateObjectiveArea()
{
    objectiveArea_ = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if
        (
            (
                p.type() == "outlet"
             || p.patch().physicalType() == "outlet"
             || p.type() == "cyclicAMI"
            )
         && (
                objectivePatch_[patchI]
            )
        )
        {
            objectiveArea_ += sum(p.magSf());
        }
    }

    reduce(objectiveArea_, sumOp<scalar>());
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getVolumeAxialDeviation()
{
    if (!axialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());

    scalar deviation = 0;
    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells(mesh_.cellZones()[index]);
        forAll(sCells, sI)
        {
            label i = sCells[sI];
            deviation += magSqr((Vel[i]&axis_)-axialVel_.target_)*Vol[i];
        }
    }
    reduce(deviation, sumOp<scalar>());
    return deviation/objectiveVolume_;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getVolumePeripheralDeviation()
{
    if (!peripheralVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    const volVectorField& cellCtr(mesh_.C());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());
    scalar deviation = 0;
    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells(mesh_.cellZones()[index]);

        forAll(sCells, sI)
        {
            label i = sCells[sI];
            vector radius = cellCtr[i] - origin_;
            vector peripheralDir = radius ^ swirlDirection_;
            peripheralDir /= mag(peripheralDir) + VSMALL;
            deviation += magSqr((Vel[i]&peripheralDir)-peripheralVel_.target_)*Vol[i];
        }
    }
    reduce(deviation, sumOp<scalar>());
    return deviation/objectiveVolume_;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getVolumeRadialDeviation()
{
    if (!radialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    const volVectorField& cellCtr(mesh_.C());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());
    scalar deviation = 0;
    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells(mesh_.cellZones()[index]);

        forAll(sCells, sI)
        {
            label i = sCells[sI];
            vector radius = cellCtr[i] - origin_;
            point hitPoint = origin_ + (radius & axis_) * axis_;
            vector hitVec = cellCtr[i] - hitPoint;
            hitVec /= mag(hitVec) + VSMALL;
            deviation += magSqr((Vel[i]&hitVec)-radialVel_.target_)*Vol[i];
        }
    }
    reduce(deviation, sumOp<scalar>());
    deviation /= objectiveVolume_;
    return deviation;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getVolumePeripheralOverAxialDeviation()
{
    if (!peripheralOverAxialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    const volVectorField& cellCtr(mesh_.C());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());
    scalar deviation = 0;
    scalar angleSum = 0;
    forAll(zoneNames_, szI)
    {
        label index = mesh_.cellZones().findZoneID(zoneNames_[szI]);
        const labelList& sCells(mesh_.cellZones()[index]);

        forAll(sCells, sI)
        {
            label i = sCells[sI];
            vector radius = cellCtr[i] - origin_;
            vector peripheralDir = radius ^ swirlDirection_;
            peripheralDir /= mag(peripheralDir) + VSMALL;
            scalar peripheralV = Vel[i]&peripheralDir;
            scalar axialV = Vel[i]&axis_;
            angleSum += mag(atan(peripheralV/(axialV + VSMALL))*Vol[i]);
        }
    }
    reduce(angleSum, sumOp<scalar>());
    deviation =
        mag(angleSum - peripheralOverAxialVel_.target_*objectiveVolume_)/objectiveVolume_;
    return deviation;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getSurfaceAxialDeviation()
{
    if (!axialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    scalar deviation = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if (isAnActivePatch(patchI))
        {
            forAll(p, fI)
            {
                vector faceVel = Vel.boundaryField()[patchI][fI];
                scalar sF = p.magSf()[fI];
                deviation += magSqr((faceVel&axis_)-axialVel_.target_)*sF;
            }
        }
    }
    reduce(deviation, sumOp<scalar>());
    return deviation/objectiveArea_;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getSurfacePeripheralDeviation()
{
    if (!peripheralVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    scalar deviation = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if (isAnActivePatch(patchI))
        {
            forAll(p, fI)
            {
                vector faceVel = Vel.boundaryField()[patchI][fI];
                vector radius = p.Cf()[fI] - origin_;
                vector peripheralDir = radius ^ swirlDirection_;
                peripheralDir /= mag(peripheralDir) + VSMALL;
                scalar sF = p.magSf()[fI];
                deviation += magSqr(mag(faceVel&peripheralDir)-peripheralVel_.target_)*sF;
            }
        }
    }
    reduce(deviation, sumOp<scalar>());

    return deviation/objectiveArea_;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getSurfaceRadialDeviation()
{
    if (!radialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    scalar deviation = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if (isAnActivePatch(patchI))
        {
            forAll(p, fI)
            {
                vector faceVel = Vel.boundaryField()[patchI][fI];
                vector radius = p.Cf()[fI] - origin_;
                scalar sF = p.magSf()[fI];
                vector peripheralDir = radius ^ swirlDirection_;
                peripheralDir /= mag(peripheralDir) + VSMALL;
                scalar peripheralV = faceVel&peripheralDir;
                scalar normalV = faceVel&p.nf()()[fI];
                normalV *= sign(p.nf()()[fI]&axis_);
                scalar radialVel =
                    mag(faceVel - peripheralV*peripheralDir - normalV*p.nf()()[fI]);
                deviation += magSqr(radialVel-radialVel_.target_)*sF;
            }
        }
    }
    reduce(deviation, sumOp<scalar>());
    deviation /= objectiveArea_;
    return deviation;
}

Foam::scalar Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
getSurfacePeripheralOverAxialDeviation()
{
    if (!peripheralOverAxialVel_.active_)
    {
        return -1;
    }
    const volVectorField& Vel(U());
    scalar deviation = 0;
    scalar angleSum = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if (isAnActivePatch(patchI))
        {
            forAll(p, fI)
            {
                vector faceVel = Vel.boundaryField()[patchI][fI];
                vector radius = p.Cf()[fI] - origin_;
                vector peripheralDir = radius ^ axis_;
                peripheralDir /= mag(peripheralDir) + VSMALL;
                scalar peripheralV = faceVel&peripheralDir;
                scalar normalV = faceVel&p.nf()()[fI];
                //For internal patches the normal could be opposite of the actual
                // flow direction
                normalV *= sign(p.nf()()[fI]&axis_);
                scalar sF = p.magSf()[fI];
                angleSum += (atan(peripheralV/(normalV + VSMALL))*sF);
            }
        }
    }
    reduce(angleSum, sumOp<scalar>());
    deviation =
        mag(angleSum - peripheralOverAxialVel_.target_*objectiveArea_)/objectiveArea_;
    return deviation;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
velocityComponentsObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    axis_(objectiveDict.lookup<vector>("axis")),
    origin_(objectiveDict.lookup<vector>("origin")),
    swirlDirection_(objectiveDict.lookupOrDefault<vector>("swirlDirection", axis_)),
    volumeBased_(objectiveDict.lookupOrDefault<bool>("volumeBased", false)),
    objectiveVolume_(0),
    objectiveArea_(0),
    axialVel_(),
    peripheralVel_(),
    radialVel_(),
    peripheralOverAxialVel_()
{
    createFiles(useAdjointFileFormat);
    if (mag(axis_) != 0 && mag(swirlDirection_) != 0)
    {
        axis_ /= mag(axis_);
        swirlDirection_ /= mag(swirlDirection_);
    }
    else
    {
        FatalErrorInFunction<< "Zero-vector axis is defined "<<exit(FatalError);
    }
    setVariables(objectiveDict);
    calculateObjectiveVolume();
    calculateObjectiveArea();
}

Foam::functionObjects::velocityComponentsObjectiveFunctionObject::
velocityComponentsObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
velocityComponentsObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::functionObjects::velocityComponentsObjectiveFunctionObject::~velocityComponentsObjectiveFunctionObject()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::velocityComponentsObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}

bool Foam::functionObjects::velocityComponentsObjectiveFunctionObject::execute()
{
    calculateObjectiveVolume();
    calculateObjectiveArea();
    if (volumeBased_)
    {
        axialVel_.deviation_ = getVolumeAxialDeviation();
        peripheralVel_.deviation_ = getVolumePeripheralDeviation();
        radialVel_.deviation_ = getVolumeRadialDeviation();
        peripheralOverAxialVel_.deviation_ = getVolumePeripheralOverAxialDeviation();
    }
    else
    {
        axialVel_.deviation_ = getSurfaceAxialDeviation();
        peripheralVel_.deviation_ = getSurfacePeripheralDeviation();
        radialVel_.deviation_ = getSurfaceRadialDeviation();
        peripheralOverAxialVel_.deviation_ = getSurfacePeripheralOverAxialDeviation();
    }

    objectiveValue_ =
        axialVel_.deviation_*axialVel_.weights_
      + radialVel_.deviation_*radialVel_.weights_
      + peripheralVel_.deviation_*peripheralVel_.weights_
      + peripheralOverAxialVel_.deviation_*peripheralOverAxialVel_.weights_;

    Info<< type() << " " << name() << " execute:" << nl
        << "Total objective value = " << objectiveValue_ << nl << endl;

    return true;
}

bool
Foam::functionObjects::velocityComponentsObjectiveFunctionObject::write()
{
    if (writeToFile() && Pstream::master())
    {
        writeTime(objFilePtr_());
        objFilePtr_()
            << tab << axialVel_.deviation_
            << tab << radialVel_.deviation_
            << tab << peripheralVel_.deviation_
            << tab << peripheralOverAxialVel_.deviation_
            << tab << objectiveValue_
            << endl;
    }

    return true;
}