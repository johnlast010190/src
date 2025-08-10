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

#include "massflowSplit/massflowSplitObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massflowSplitObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        massflowSplitObjectiveFunctionObject,
        dictionary
    );
}
}
// * * * * * * * * * * * * * *Private Member Functions * * * * * * * * * * * //
void Foam::functionObjects::massflowSplitObjectiveFunctionObject::writeFileHeader
(
    Ostream& os
)
{
    writeCommented(os, "Time actualMassFlowSplit ( ");

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            if (Pstream::master())
            {
                writeDelimited(os, mesh_.boundary()[patchI].name());
            }
        }
    }
    writeDelimited(os, ") objectiveValue");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massflowSplitObjectiveFunctionObject::
massflowSplitObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    massflowFraction_(mesh_.boundary().size(), 0.0)
{
    read(objectiveDict);
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::massflowSplitObjectiveFunctionObject::
massflowSplitObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    massflowSplitObjectiveFunctionObject(name, runTime, objectiveDict, false)
{
    read(objectiveDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massflowSplitObjectiveFunctionObject::
~massflowSplitObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::massflowSplitObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    objectivePatch_ = false;

    if (dict.found("flowRateFractions"))
    {
        List<Tuple2<word, scalar>> patchFracs
        (
            dict.lookup("flowRateFractions")
        );

        forAll(patchFracs, pfI)
        {
            label patchI
                = mesh_.boundaryMesh().findPatchID(patchFracs[pfI].first());

            if (patchI != -1)
            {
                objectivePatch_[patchI] = true;
                massflowFraction_[patchI] = patchFracs[pfI].second();
            }
            else
            {
                FatalErrorInFunction
                    << "Patch " << patchFracs[pfI].first()
                    << " specified by flowRateFractions not found."
                    << exit(FatalError);
            }
        }
    }
    else
    {
        //equal fractions for each outlet
        scalarList massflowFraction_(mesh_.boundary().size(), 0.0);

        forAll(mesh_.boundary(), patchI)
        {
            const fvPatch& p = mesh_.boundary()[patchI];

            if (p.patch().physicalType() == "outlet" || p.type() == "outlet")
            {
                massflowFraction_[patchI] = 1;
                objectivePatch_[patchI] = true;
            }
        }

    }

    scalar fracSum(sum(massflowFraction_));

    if (mag(1.0 - fracSum) > SMALL)
    {
        WarningInFunction
            << "Sum of mass flow target fractions do not add to one."
            << endl;

    }

    return true;
}


bool
Foam::functionObjects::massflowSplitObjectiveFunctionObject::execute()
{
    scalar totalInletFlowRate = inletFlowRate();

    objectiveValue_ = 0;

    scalar fracSum(sum(massflowFraction_));

    if (mag(1.0 - fracSum) > SMALL)
    {
        WarningInFunction
            << "Sum of mass flow target fractions do not add to one."
            << endl;
    }

    if (writeToFile() && Pstream::master())
    {
        writeTime(objFilePtr_());
    }

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            scalar patchFlow = gSum(phi().boundaryField()[patchI]);

            if (writeToFile() && Pstream::master())
            {
                writeDelimited(objFilePtr_(), mag(patchFlow/totalInletFlowRate));
            }

            scalar dMk;
            if (patchFlow > 0)
            {
                dMk = patchFlow
                    - massflowFraction_[patchI]*totalInletFlowRate;
            }
            else
            {
                dMk = -patchFlow
                    + massflowFraction_[patchI]*totalInletFlowRate;
            }
            objectiveValue_ += 0.5*sqr(dMk);
        }
    }

    Info<< type() << " " << name() << " execute:" << nl
        << "Target mass flow split RMS = " << objectiveValue_
        << " [kg/s]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::massflowSplitObjectiveFunctionObject::write()
{
    if (writeToFile() && Pstream::master())
    {
            objFilePtr_()
                << tab << objectiveValue_
                << endl;
    }
    return true;
}


// ************************************************************************* //
