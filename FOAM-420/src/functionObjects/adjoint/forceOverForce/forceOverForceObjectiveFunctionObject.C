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

#include "forceOverForce/forceOverForceObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceOverForceObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        forceOverForceObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::vector
Foam::functionObjects::forceOverForceObjectiveFunctionObject::
totalForce() const
{
    tmp<volScalarField> P = objectiveFunctionObject::P();

    const volScalarField::Boundary& Pw
        = P().boundaryField();

    tmp<volSymmTensorField> tau = devRhoReff();

    const volSymmTensorField::Boundary& tauw
        = tau().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    vector totalForce = vector::zero;

    forAll(objectivePatch_, pI)
    {
        if (objectivePatch_[pI] && isA<wallFvPatch>(mesh_.boundary()[pI]))
        {
            vectorField pf( Sfb[pI]*Pw[pI] );
            vectorField vf( Sfb[pI] & tauw[pI] );

            totalForce += sum(pf + vf);
        }
    }

    reduce(totalForce, sumOp<vector>());

    return totalForce;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceOverForceObjectiveFunctionObject::
forceOverForceObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    numForceDirection_(objectiveDict.lookup("numeratorForceDirection")),
    denForceDirection_(objectiveDict.lookup("denominatorForceDirection"))
{
    createFiles(useAdjointFileFormat);

    numForceDirection_ /= mag(numForceDirection_);
    denForceDirection_ /= mag(denForceDirection_);
}


Foam::functionObjects::forceOverForceObjectiveFunctionObject::
forceOverForceObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    forceOverForceObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceOverForceObjectiveFunctionObject::
~forceOverForceObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::forceOverForceObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    numForceDirection_ = vector(dict.lookup("numeratorForceDirection"));
    denForceDirection_ = vector(dict.lookup("denominatorForceDirection"));

    return true;
}


bool
Foam::functionObjects::forceOverForceObjectiveFunctionObject::
execute()
{
    vector totalF = totalForce();

    scalar numForce = totalF & numForceDirection_;
    scalar denForce = totalF & denForceDirection_;

    objectiveValue_ = numForce/denForce;

    Info<< type() << " " << name() << " execute:" << nl
        << "Force over force = " << objectiveValue_ << " []" << nl << endl;

    return true;
}


bool
Foam::functionObjects::forceOverForceObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
