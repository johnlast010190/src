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

#include "composite/compositeObjectiveFunctionObject.H"
//#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(compositeObjectiveFunctionObject, 0);
    //addToRunTimeSelectionTable
    //(
    //    functionObject,
    //    compositeObjectiveFunctionObject,
    //    dictionary
    //);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::compositeObjectiveFunctionObject::
compositeObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    compositeObjectives_(),
    normalisedCompositeWeights_(0)
{
    createFiles(useAdjointFileFormat);

    dictionary& refDict = const_cast<dictionary&>(objectiveDict);

    label nObjectives = refDict.subDict("objectives").size();

    compositeObjectives_.setSize(nObjectives);
    normalisedCompositeWeights_.setSize(nObjectives);

    label i = 0;
    forAllConstIter(dictionary, refDict.subDict("objectives"), iter)
    {
        compositeObjectives_.set
        (
            i++,
            functionObject::New
            (
                iter().dict().dictName(),
                runTime,
                iter().dict()
            )
        );
    }
}


Foam::functionObjects::compositeObjectiveFunctionObject::
compositeObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    compositeObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::compositeObjectiveFunctionObject::
~compositeObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::compositeObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::compositeObjectiveFunctionObject::
execute()
{
    objectiveValue_ = 0.0;

    forAll(compositeObjectives_, oI)
    {
        compositeObjectives_[oI].execute();

        objectiveFunctionObject& objFoI =
            refCast<objectiveFunctionObject>
            (
                compositeObjectives_[oI]
            );
        objectiveValue_ +=
            objFoI.objectiveValue()*normalisedCompositeWeights_[oI];
    }

    Info<< type() << " " << name() << " execute:" << nl
        << "Combined value = " << objectiveValue_ << " []" << nl << endl;

    return true;
}


bool
Foam::functionObjects::compositeObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
