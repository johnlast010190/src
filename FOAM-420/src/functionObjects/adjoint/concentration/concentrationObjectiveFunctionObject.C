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

#include "concentration/concentrationObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(concentrationObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        concentrationObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * *Private Member Functions * * * * * * * * * * * //

void
Foam::functionObjects::concentrationObjectiveFunctionObject::
readData(const dictionary& dict)
{
    List<Tuple2<word, scalar>> patchTarget
    (
        dict.lookup("psiTarget")
    );

    forAll(patchTarget, pfI)
    {
        label patchI
            = mesh_.boundaryMesh().findPatchID(patchTarget[pfI].first());
        if (patchI != -1)
        {
            objectivePatch_[patchI] = true;
            psiTar_[patchI] = patchTarget[pfI].second();
        }
        else
        {
            FatalErrorInFunction
                << "Patch " << patchTarget[pfI].first()
                << " specified by psiTarget not found."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::concentrationObjectiveFunctionObject::
concentrationObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    psiTar_(mesh_.boundary().size(), 0.0)
{
    createFiles(useAdjointFileFormat);

    readData(objectiveDict);
}


Foam::functionObjects::concentrationObjectiveFunctionObject::
concentrationObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    concentrationObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::concentrationObjectiveFunctionObject::
~concentrationObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::concentrationObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    readData(dict);

    return true;
}


bool
Foam::functionObjects::concentrationObjectiveFunctionObject::
execute()
{
    scalar cost = 0;
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];
        if
        (
            (p.type() == "outlet" || p.patch().physicalType() == "outlet")
            && (objectivePatch_[patchI])
        )
        {
            cost += sum
            (
                0.5*p.magSf()
               *mag
                (
                    psi().boundaryField()[patchI]
                  - psiTar_[patchI]
                )
            );
        }
    }

    reduce(cost, sumOp<scalar>());

    objectiveValue_ = cost;

    Info<< type() << " " << name() << " execute:" << nl
        << "Deviation from target concentration = " << objectiveValue_
        << " []" << nl << endl;

    return true;
}


bool
Foam::functionObjects::concentrationObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
