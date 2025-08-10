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
    (c) 2017 OpenCFD Ltd.

Description
    Basic tests of IOobjectList
\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/timeSelector.H"
#include "db/IOobjectList/IOobjectList.H"
#include "primitives/strings/lists/hashedWordList.H"
#include "fields/volFields/volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("re", "wordReList");

    // timeSelector::addOptions();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"

    wordReList matcher;
    if (args.optionFound("re"))
    {
        matcher = args.optionReadList<wordRe>("re");
        Info<<"limit names: " << matcher << nl;

    }

    const hashedWordList subsetTypes
    {
        volScalarField::typeName,
        volScalarField::Internal::typeName,
        volVectorField::typeName,
    };


    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        // Objects at this time
        IOobjectList objects(runTime, runTime.timeName());
        HashTable<wordHashSet> classes =
        (
            matcher.size()
          ? objects.classes(matcher)
          : objects.classes()
        );

        Info<< "Time: " << runTime.timeName() << nl;

        Info<<"Name:    " << flatOutput(objects.sortedNames()) << nl
            <<"Objects: " << objects << nl
            <<"Classes: " << classes << nl;

        classes.filterKeys(subsetTypes);
        Info<<"only retain: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        classes = objects.classes();
        classes.erase(subsetTypes);
        Info<<"remove: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
