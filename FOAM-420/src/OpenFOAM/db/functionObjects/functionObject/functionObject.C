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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/functionObjects/functionObject/functionObject.H"
#include "db/dictionary/dictionary.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(functionObject, "functionObject", 0);
    defineRunTimeSelectionTable(functionObject, dictionary);
}

bool Foam::functionObject::postProcess(false);

Foam::word Foam::functionObject::outputPrefix("postProcessing");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObject::scopedName(const word& name) const
{
    return name_ + ":" + name;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObject::functionObject(const word& name)
:
    name_(name),
    log(postProcess)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObject> Foam::functionObject::New
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
{
    const word functionType(dict.lookup("type"));

    if (debug)
    {
        Info<< "Selecting function " << functionType << endl;
    }

    if (dict.found("functionObjectLibs"))
    {
        const_cast<Time&>(runTime).libs().open
        (
            dict,
            "functionObjectLibs",
            dictionaryConstructorTable_()
        );
    }
    else
    {
        const_cast<Time&>(runTime).libs().open
        (
            dict,
            "libs",
            dictionaryConstructorTable_()
        );
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTable_().find(functionType);


    if (cstrIter == dictionaryConstructorTable_().end())
    {
        Switch exitIfNotFound
        (
            dict.lookupOrDefault<Switch>("exitIfNotFound", true)
        );

        if (exitIfNotFound)
        {
            FatalErrorInFunction
                << "Unknown function type "
                << functionType << nl << nl
                << exit(FatalError);
        }
        else
        {
            WarningInFunction
                << "Unknown function type "
                << functionType << nl << nl
                << "Continuing with null function." << endl
                << endl;

            cstrIter = dictionaryConstructorTable_().find("null");
        }
    }

    return autoPtr<functionObject>(cstrIter->second(name, runTime, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::functionObject::name() const
{
    return name_;
}


bool Foam::functionObject::read(const dictionary& dict)
{
    if (!postProcess)
    {
        log = dict.lookupOrDefault<Switch>("log", true);
    }

    return true;
}


bool Foam::functionObject::execute(const label)
{
    return true;
}


bool Foam::functionObject::end()
{
    return true;
}


bool Foam::functionObject::adjustTimeStep()
{
    return false;
}


bool Foam::functionObject::filesModified() const
{
    return false;
}


void Foam::functionObject::updateMesh(const mapPolyMesh&)
{}


void Foam::functionObject::movePoints(const polyMesh&)
{}


// ************************************************************************* //
