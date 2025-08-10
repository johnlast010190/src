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

Application
    foamList

Description
    Print the table of contents of selectable switches, classes etc. in the
    OpenFOAM libraries

Usage
    \b foamList [OPTION]

    Options:
      - \par -switches
        Print the DebugSwitches, InfoSwitches and OptimisationSwitches

      - \par -registeredSwitches
        Print the registered DebugSwitches, InfoSwitches and
        OptimisationSwitches supporting run-time modification

      - \par -unset
        print switches declared in libraries but not set in etc/controlDict

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/dictionary/dictionary.H"
#include "global/debug/simpleObjectRegistry.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOobject/IOobject.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "global/etcFiles/etcFiles.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "db/functionObjects/functionObject/functionObject.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

using namespace Foam;

void listSwitches
(
    const wordList& debugSwitches,
    const wordList& infoSwitches,
    const wordList& optSwitches,
    const bool unset
)
{
    if (unset)
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        dictionary controlDict;
        forAllReverse(controlDictFiles, cdfi)
        {
            controlDict.merge(dictionary(IFstream(controlDictFiles[cdfi])()));
        }

        wordHashSet controlDictDebug
        (
            controlDict.subDict("DebugSwitches").sortedToc()
        );

        wordHashSet controlDictInfo
        (
            controlDict.subDict("InfoSwitches").sortedToc()
        );

        wordHashSet controlDictOpt
        (
            controlDict.subDict("OptimisationSwitches").sortedToc()
        );


        IOobject::writeDivider(Info);

        wordHashSet hashset;
        hashset = debugSwitches;
        hashset -= controlDictDebug;
        Info<< "Unset DebugSwitches" << hashset.sortedToc() << endl;

        hashset = infoSwitches;
        hashset -= controlDictInfo;
        Info<< "Unset InfoSwitches" << hashset.sortedToc() << endl;

        hashset = optSwitches;
        hashset -= controlDictOpt;
        Info<< "Unset OptimisationSwitches" << hashset.sortedToc() << endl;
    }
    else
    {
        IOobject::writeDivider(Info);
        Info<< "DebugSwitches" << debugSwitches << endl;
        Info<< "InfoSwitches" << infoSwitches << endl;
        Info<< "OptimisationSwitches" << optSwitches << endl;
    }
}


void listSwitches(const argList& args)
{
    if (args.optionFound("registeredSwitches"))
    {
        listSwitches
        (
            debug::debugObjects().sortedToc(),
            debug::infoObjects().sortedToc(),
            debug::optimisationObjects().sortedToc(),
            args.optionFound("unset")
        );
    }
    else
    {
        listSwitches
        (
            debug::debugSwitches().sortedToc(),
            debug::infoSwitches().sortedToc(),
            debug::optimisationSwitches().sortedToc(),
            args.optionFound("unset")
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "switches",
        "List switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "registeredSwitches",
        "List switches registered for run-time modification"
    );
    argList::addBoolOption
    (
        "unset",
        "List switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "scalarBCs",
        "List scalar field boundary conditions (fvPatchField<scalar>)"
    );
    argList::addBoolOption
    (
        "vectorBCs",
        "List vector field boundary conditions (fvPatchField<vector>)"
    );
    argList::addBoolOption
    (
        "functionObjects",
        "List functionObjects"
    );
    argList::addBoolOption
    (
        "fvOptions",
        "List fvOptions"
    );
    argList::addBoolOption
    (
        "incompressibleTurbulenceModels",
        "List incompressible turbulenceModels"
    );
    argList::addBoolOption
    (
        "compressibleTurbulenceModels",
        "List compressible turbulenceModels"
    );

    argList args(argc, argv);

    if (!args.options().size())
    {
        args.printUsage();
        return 0;
    }

    if
    (
        args.optionFound("switches")
     || args.optionFound("registeredSwitches")
    )
    {
        listSwitches(args);
    }

    if (args.optionFound("scalarBCs"))
    {
        Info<< "scalarBCs" << nl;
        printCtorTableKeys
        (
            Info, fvPatchField<scalar>::dictionaryConstructorTable_()
        );
        Info<< endl;
    }
    if (args.optionFound("vectorBCs"))
    {
        Info<< "vectorBCs" << nl;
        printCtorTableKeys
        (
            Info, fvPatchField<vector>::dictionaryConstructorTable_()
        );
        Info<< endl;
    }

    if (args.optionFound("functionObjects"))
    {
        Info<< "functionObjects" << nl;
        printCtorTableKeys
        (
            Info, functionObject::dictionaryConstructorTable_()
        );
        Info<< endl;
    }

    if (args.optionFound("fvOptions"))
    {
        Info<< "fvOptions" << nl;
        printCtorTableKeys
        (
            Info, fv::option::dictionaryConstructorTable_()
        );
        Info<< endl;
    }

    if (args.optionFound("incompressibleTurbulenceModels"))
    {
        Info<< "Turbulence models" << nl;
        printCtorTableKeys
        (
            Info, incompressible::turbulenceModel::dictionaryConstructorTable_()
        );
        Info<< endl;

        Info<< "RAS models" << nl;
        printCtorTableKeys
        (
            Info, incompressible::RASModel::dictionaryConstructorTable_()
        );
        Info<< endl;

        Info<< "LES models" << nl;
        printCtorTableKeys
        (
            Info, incompressible::LESModel::dictionaryConstructorTable_()
        );
        Info<< endl;
    }

    if (args.optionFound("compressibleTurbulenceModels"))
    {
        Info<< "Turbulence models" << nl;
        printCtorTableKeys
        (
            Info, compressible::turbulenceModel::dictionaryConstructorTable_()
        );
        Info<< endl;

        Info<< "RAS models" << nl;
        printCtorTableKeys
        (
            Info, compressible::RASModel::dictionaryConstructorTable_()
        );
        Info<< endl;

        Info<< "LES models" << nl;
        printCtorTableKeys
        (
            Info, compressible::LESModel::dictionaryConstructorTable_()
        );
        Info<< endl;
    }


    return 0;
}


// ************************************************************************* //
