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

#include "torque/torqueObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(torqueObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        torqueObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::functionObjects::torqueObjectiveFunctionObject::
calcTorque() const
{
    tmp<volScalarField> P = objectiveFunctionObject::P();

    const volScalarField::Boundary& Pw
        = P().boundaryField();

    tmp<volSymmTensorField> tau = devRhoReff();

    const volSymmTensorField::Boundary& tauw
        = tau().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    vector totalTorque = vector::zero;

    forAll(objectivePatch_, pI)
    {
        if (objectivePatch_[pI] && isA<wallFvPatch>(mesh_.boundary()[pI]))
        {
            vectorField pf( Sfb[pI]*Pw[pI] );
            vectorField vf( Sfb[pI] & tauw[pI] );
            vectorField ff( pf + vf );

            vectorField CfpI = mesh_.C().boundaryField()[pI];
            pointField hitPoint
            (
                cAxisRotation_
                + ((CfpI - cAxisRotation_) & torqueDirection_) * torqueDirection_
            );

            CfpI -= hitPoint;

            vectorField torque( CfpI^(ff) );

            totalTorque += sum(torque);
        }
    }

    reduce(totalTorque, sumOp<vector>());

    return totalTorque;
}


Foam::tmp<vectorField>
Foam::functionObjects::torqueObjectiveFunctionObject::
torqueWall
(
    const fvPatchVectorField& pf
) const
{
    tmp<vectorField> tvMod(new vectorField(pf.size(), vector::zero));
    vectorField& vMod = tvMod.ref();

    label index = pf.patch().index();

    forAll(mesh_.boundary()[index], pfI)
    {
        vector CfI = mesh_.C().boundaryField()[index][pfI];
        vector rCfI = CfI - cAxisRotation_;

        vMod[pfI].x() = torqueDirection_.y() * rCfI.z() -
                        torqueDirection_.z() * rCfI.y();

        vMod[pfI].y() = torqueDirection_.z() * rCfI.x() -
                        torqueDirection_.x() * rCfI.z();

        vMod[pfI].z() = torqueDirection_.x() * rCfI.y() -
                        torqueDirection_.y() * rCfI.x();
    }

    vMod *= sign(-objectiveValue_);

    return tvMod;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::torqueObjectiveFunctionObject::
torqueObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    torqueDirection_(objectiveDict.lookup("torqueDirection")),
    cAxisRotation_(objectiveDict.lookup("centreAxisRotation"))
{
    createFiles(useAdjointFileFormat);

    torqueDirection_ /= mag(torqueDirection_);

    vector totalTorque = calcTorque();
    objectiveValue_ = totalTorque & torqueDirection_;
}


Foam::functionObjects::torqueObjectiveFunctionObject::
torqueObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    torqueObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::torqueObjectiveFunctionObject::
~torqueObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::torqueObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    torqueDirection_ = vector(dict.lookup("torqueDirection"));
    torqueDirection_ /= mag(torqueDirection_);
    cAxisRotation_ = vector(dict.lookup("centreAxisRotation"));

    return true;
}


bool
Foam::functionObjects::torqueObjectiveFunctionObject::
execute()
{
    vector totalTorque = calcTorque();

    objectiveValue_ = (totalTorque & torqueDirection_);

    Info<< type() << " " << name() << " execute:" << nl
        << "Torque = " << objectiveValue_ << " [Nm]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::torqueObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
