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

#include "temperature/temperatureObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(temperatureObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        temperatureObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::temperatureObjectiveFunctionObject::
temperatureObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict))
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::temperatureObjectiveFunctionObject::
temperatureObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    temperatureObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::temperatureObjectiveFunctionObject::
~temperatureObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::temperatureObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::temperatureObjectiveFunctionObject::
execute()
{
    const volScalarField& T = temperature();

    const volScalarField::Boundary& Tp = T.boundaryField();

    const surfaceScalarField::Boundary& Afp =
        mesh_.magSf().boundaryField();

    scalar totalT = 0.0;
    scalar surface = 0.0;

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            totalT += sum(Tp[patchI]*Afp[patchI]);
            surface += sum(Afp[patchI]);
        }
    }
    reduce(totalT, sumOp<scalar>());
    reduce(surface, sumOp<scalar>());
    objectiveValue_ = totalT/surface;

    Info<< type() << " " << name() << " execute:" << nl
        << "Mean temperature = " << objectiveValue_
        << " [C]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::temperatureObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
