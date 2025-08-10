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

#include "swirl/swirlObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(swirlObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        swirlObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * *Protected Member Functions * * * * * * * * * * * //

Foam::scalar
Foam::functionObjects::swirlObjectiveFunctionObject::calcSwirl() const
{
    scalar swirl = 0;

    tmp<volScalarField> rho = this->rho();
    const volVectorField& Vel(U());
    const volVectorField& cC(mesh_.C());
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());

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
            vector c = cC[i] - axisPoint_;

            point hitPoint = axisPoint_
                + (c & swirlDirection_) * swirlDirection_;
            point hitVec = cC[i] - hitPoint;
            scalar dist = mag(hitVec);
            point tangentialDir = swirlDirection_ ^ (hitVec / (dist + SMALL));

            swirl
                += Vol[i]*rho()[i]*
                    (Vel[i]&swirlDirection_)*
                    ((hitVec^((Vel[i]&tangentialDir)*tangentialDir))&swirlDirection_);
        }
    }

    reduce(swirl, sumOp<scalar>());

    return swirl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::swirlObjectiveFunctionObject::
swirlObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    axisPoint_
    (
        objectiveDict.lookupOrDefault<vector>
        (
            "axisPoint",
            vector(0.0, 0.0, 0.0)
        )
    ),
    swirlDirection_
    (
        objectiveDict.lookupOrDefault<vector>
        (
            "swirlDirection",
            vector(0.0, 1.0, 0.0)
        )
    )
{
    createFiles(useAdjointFileFormat);
    swirlDirection_ /= mag(swirlDirection_);
}


Foam::functionObjects::swirlObjectiveFunctionObject::
swirlObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    swirlObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::swirlObjectiveFunctionObject::
~swirlObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::swirlObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::swirlObjectiveFunctionObject::execute()
{
    scalar mass = 0;

    tmp<volScalarField> rho = this->rho();
    const DimensionedField<scalar, volMesh>& Vol(mesh_.V());

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
            mass += Vol[i]*rho()[i];
        }
    }

    objectiveValue_ = mag(calcSwirl());
    reduce(mass, sumOp<scalar>());
    objectiveValue_ /= mass;

    Info<< type() << " " << name() << " execute:" << nl
        << "Mean swirl = " << objectiveValue_ << " [m3/s2]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::swirlObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
