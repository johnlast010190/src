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

#include "force/forceObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        forceObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceObjectiveFunctionObject::
forceObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    forceDirection_(objectiveDict.lookup("forceDirection"))
{
    createFiles(useAdjointFileFormat);
    forceDirection_ /= mag(forceDirection_);
}


Foam::functionObjects::forceObjectiveFunctionObject::
forceObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    forceObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceObjectiveFunctionObject::
~forceObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::forceObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::forceObjectiveFunctionObject::execute()
{
    tmp<volScalarField> pp = objectiveFunctionObject::P();
    const volScalarField::Boundary& Pw
        = pp().boundaryField();

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
    objectiveValue_ = (totalForce & forceDirection_);

    Info<< type() << " " << name() << " execute:" << nl
        << "Total force = " << objectiveValue_ << " [N]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::forceObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
