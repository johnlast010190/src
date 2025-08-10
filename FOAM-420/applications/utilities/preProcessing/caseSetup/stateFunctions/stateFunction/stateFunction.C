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
    (c) 2016-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "meshes/polyMesh/polyMesh.H"
#include "global/etcFiles/etcFiles.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/error/error.H"
#include "stateFunctions/global/globalState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(stateFunction, 0);
defineRunTimeSelectionTable(stateFunction, dictionary);

template<>
const char* NamedEnum<stateFunction::regionType, 5>::names[] =
{
    "fluid",
    "solid",
    "surface",
    "none",
    "multipointAdjoint"
};

template<>
const char* NamedEnum<stateFunction::timeType, 3>::names[] =
{
    "steady",
    "transient",
    "none"
};

template<>
const char* NamedEnum<stateFunction::turbulenceType, 4>::names[] =
{
    "RAS",
    "LES",
    "laminar",
    "none"
};

template<>
const char* NamedEnum<stateFunction::compressibilityType, 3>::names[] =
{
    "incompressible",
    "compressible",
    "none"
};

template<>
const char* NamedEnum<stateFunction::customSolver, 4>::names[] =
{
    "standard",
    "speciesMixing",
    "tempFoam",
    "none"
};

template<>
const char* NamedEnum<stateFunction::matrixType, 2>::names[] =
{
    "segregated",
    "coupled"
};

template<>
const char* NamedEnum<stateFunction::meshType, 3>::names[] =
{
    "volume",
    "point",
    "surface"
};

template<>
const char* NamedEnum<stateFunction::fieldType, 5>::names[] =
{
    "scalar",
    "vector",
    "tensor",
    "symmTensor",
    "sphericalTensor"
};

const NamedEnum<stateFunction::regionType, 5> stateFunction::regionTypeNames_;

const NamedEnum<stateFunction::timeType, 3> stateFunction::timeTypeNames_;

const NamedEnum<stateFunction::turbulenceType, 4>
    stateFunction::turbulenceTypeNames_;

const NamedEnum<stateFunction::compressibilityType, 3>
    stateFunction::compressibilityTypeNames_;

const NamedEnum<stateFunction::customSolver, 4>
    stateFunction::customSolverNames_;

const NamedEnum<stateFunction::matrixType, 2>
    stateFunction::matrixTypeNames_;

const NamedEnum<stateFunction::meshType, 3>
    stateFunction::meshTypeNames_;

const NamedEnum<stateFunction::fieldType, 5>
    stateFunction::fieldTypeNames_;


// * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * *//

fileName stateFunction::fileDir
(
    const Time& time,
    const fileName& outputDir,
    const word& region
)
{
    //Does not take care of parallel execution effects
    //should be taken care of via ../ in outputDir
    fileName outPath(time.time().path()/outputDir);

    // Append mesh name if not default region or nullptr
    if (region != word::null && region != polyMesh::defaultRegion)
    {
        outPath = outPath/region;
    }

    outPath.expand();
    outPath.toAbsolute();

    return outPath;
}


void stateFunction::writeDirDicts
(
    const objectRegistry& obr,
    const fileName& region,
    const fileName& local,
    const dictionary& dicts,
    const wordList& excludedDicts,
    const wordList& includedDicts
)
{
    word lclreg(region);
    if (lclreg == polyMesh::defaultRegion)
    {
        lclreg = word::null;
    }

    forAllConstIter(dictionary, dicts, iter)
    {
        if
        (
            includedDicts.size()
          ? includedDicts.found(iter().keyword())
          : !excludedDicts.found(iter().keyword())
        )
        {
            IOdictionary odict
            (
                IOobject
                (
                    iter().keyword(),
                    local,
                    lclreg,
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                iter().dict()
            );
            if (lclreg == word::null)
            {
                Info<< "Writing " << word(local)/iter().keyword() << endl;
            }
            else
            {
                Info<< "Writing "
                    << word(local)/lclreg/iter().keyword() << endl;
            }
            odict.regIOobject::writeObject
            (
                obr.time().writeFormat(),
                IOstream::currentVersion,
                IOstream::UNCOMPRESSED,
                true
            );
        }
    }
}

bool stateFunction::wordListFind(const wordList& wl, const word& w)
{
    forAll(wl, i)
    {
        if (wl[i] == w)
        {
            return true;
        }
    }

    return false;
}


bool stateFunction::wordListMatch(const wordList& a, const wordList& b)
{
    // check size
    if (a.size() != b.size())
    {
        return false;
    }

    //check duplications
    forAll(a, i)
    {
        for (label j = i + 1; j < a.size() ; j++)
        {
            if (a[i] == a[j])
            {
                WarningInFunction
                     << "List contains duplicate entries, cannot be uniquely matched."
                     << endl;
                return false;
            }
        }
    }

    //compare elements - N squared
    bool listMatch = true;

    forAll(a, i)
    {
        bool matched = false;

        forAll(b, j)
        {
            if (a[i] == b[j])
            {
                matched = true;
            }
        }

        if (!matched) return false;
    }

    return listMatch;
}


Xfer<dictionary> stateFunction::etcDictionary
(
    const fileName& fn,
    bool mustfind
)
{
    dictionary etcDict;

    fileNameList etcDictFiles(findEtcFiles(fn, mustfind));

    forAllReverse(etcDictFiles, cdfi)
    {
        IFstream ifs(etcDictFiles[cdfi]);

        if (!ifs.good())
        {
            FatalIOErrorInFunction
            (
                ifs
            )   << "Cannot open " << fn;
        }
        etcDict.merge(dictionary(ifs));
    }

    return etcDict.xfer();
}


wordList stateFunction::extractModules
(
    const dictionary& modules,
    wordList& stateList
)
{
    DynamicList<word> coreState(1);
    DynamicList<word> moduleIds(0);

    forAll(stateList, i)
    {
        if (modules.found(stateList[i]))
        {
            moduleIds.append(stateList[i]);
        }
        else
        {
            coreState.append(stateList[i]);
        }
    }

    stateList = coreState;

    return wordList(moduleIds, true);
}


word stateFunction::stateIdentifier
(
    const dictionary& input,
    const dictionary& defaults
)
{
    word stateId = stateFunctions::globalState::typeName;

    if (input.found("state"))
    {
        //extract modules from state input
        wordList stateList(input.lookup("state"));
        wordList moduleIds = extractModules
        (
            defaults.subDict("solvers").subDict("modules"),
            stateList
        );

        // match state input to available state
        PtrList<entry> stateRegister(defaults.lookup("states"));
        stateId = matchState(stateRegister, stateList);
    }


    return stateId;
}


word stateFunction::matchState
(
    const PtrList<entry>& stateRegister,
    const wordList& stateList
)
{
    forAll(stateRegister, i)
    {
        if (wordListMatch(stateList, wordList(stateRegister[i].stream())))
        {
            return stateRegister[i].keyword();
        }
    }

    FatalErrorInFunction
        << "State (" <<  stateList << ") not found." << nl
        << "Available states are: " << stateRegister << exit(FatalError);

    return word::null;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void stateFunction::insertModules
(
    const dictionary& defaults,
    const wordList& moduleIds,
    word stateId,
    dictionary& targetDict
)
{
    const dictionary& moduleDicts
        = defaults.subDict("solvers").subDict("modules");

    if (moduleIds.size() > 0)
    {
        Info<< " - modules: " << endl;
    }

    forAll(moduleIds, mI)
    {
        Info<< " -- " << moduleIds[mI] << endl;

        //check requirements AND/OR/NOT
        dictionary mDict = moduleDicts.subDict(moduleIds[mI]);

        if (mDict.found("requirements"))
        {
            const dictionary& requirements = mDict.subDict("requirements");

            if (requirements.found("AND"))
            {
                const wordList AND = requirements.lookup("AND");

                forAll(AND, rI)
                {
                    if (!wordListFind(moduleIds, AND[rI]) && stateId != AND[rI])
                    {
                        FatalError << "Module dependencies: Requisite module "
                             << AND[rI]
                             << " not found. Module " << moduleIds[mI]
                             << " must be used in conjunction with ALL"
                             << " of the following modules: " << nl
                             << AND << endl;
                        break;
                    }
                }
            }
            if (requirements.found("OR"))
            {
                const wordList OR = requirements.lookup("OR");
                DynamicList<word> modulesFoundInOr(0);

                bool found = false;
                forAll(OR, rI)
                {
                    if (wordListFind(moduleIds, OR[rI]) || stateId == OR[rI])
                    {
                        found = true;
                        modulesFoundInOr.append(OR[rI]);
                    }
                }

                if (!found)
                {
                    FatalError << "Module dependencies: One or more modules "
                         << "from the following list must be used in "
                         << "conjunction with " << moduleIds[mI] << " module: "
                         << OR << exit(FatalError);

                }

                //conditionals
                // add entries based on OR modules
                if (requirements.found("conditional"))
                {
                    modulesFoundInOr.shrink();
                    forAll(modulesFoundInOr, rI)
                    {
                        if
                        (
                            requirements.subDict("conditional")
                            .found(modulesFoundInOr[rI])
                        )
                        {
                            dictionary conditional;
                            conditional = requirements.subDict("conditional")
                                .subDict(modulesFoundInOr[rI]);
                            mDict.merge(conditional);
                        }
                    }
                }

            }
            if (requirements.found("NOT"))
            {
                const wordList NOT = requirements.lookup("NOT");

                forAll(NOT, rI)
                {
                    if (wordListFind(moduleIds, NOT[rI]) || stateId == NOT[rI])
                    {
                        FatalError << "Module dependencies: Prohibited module "
                             << NOT[rI]
                             << " found. Module " << moduleIds[mI]
                             << " cannot be used in conjunction with ANY"
                             << " of the following modules: " << nl
                             << NOT << endl;
                        break;
                    }
                }
            }
            mDict.remove("requirements");
        }


        //add module dicts to stateDict
        targetDict.merge(mDict);
    }

    Info<< endl;

}


void stateFunction::insertProgramaticSettings
(
    const dictionary& defaults,
    const dictionary& input
)
{
    // Insert buoyant settins if they are available
    bool isBuoyant = false;
    bool isMultiphase = false;
    const wordList materials = input.lookup<wordList>("materials");
    const wordList states = input.lookup<wordList>("state");
    if (input.isDict("materialProperties"))
    {
        const dictionary& matInputDict = input.subDict("materialProperties");
        if (materials.size() > 1)
        {
            isMultiphase =
                (
                    matInputDict.lookupOrDefault<word>("materialType", "fluid")
                 == "multiphase"
                );
            if (matInputDict.found("buoyancyModel"))
            {
                isBuoyant = true;
            }
            else if (isMultiphase)
            {
                isBuoyant = true;
            }
            else if (matInputDict.isDict(materials.first()) && matInputDict.subDict(materials.first()).found("buoyancyModel"))
            {
                // The model has to be specified for all species
                // or on the phase level
                isBuoyant = true;
                forAll(materials, speciesI)
                {
                    if (!matInputDict.subDict(materials[speciesI]).found("buoyancyModel"))
                    {
                        isBuoyant = false;
                    }
                }
            }
        }
        else if (materials.size() == 1)
        {
            if (matInputDict.isDict(materials.first()))
            {
                if (matInputDict.subDict(materials.first()).found("buoyancyModel"))
                {
                    isBuoyant = true;
                }
            }
        }
        if (!materials.size())
        {
            FatalErrorInFunction
                << "At least one material has to exist"
                << exit(FatalError);
        }
    }
    if (states.found("buoyant"))
    {
        isBuoyant = true;
    }
    if (isBuoyant && defaults.isDict("buoyantSteady"))
    {
        if (states.found("steady"))
        {
            const_cast<dictionary&>(defaults).merge
            (
                defaults.subDict("buoyantSteady")
            );
        }
        else if (states.found("transient"))
        {
            const_cast<dictionary&>(defaults).merge
            (
                defaults.subDict("buoyantTransient")
            );
        }
    }

    // Insert multicomponent settings
    bool isMulticomponent = false;
    if (isMultiphase)
    {
        forAll(materials, mati)
        {
            if
            (
                input.subOrEmptyDict
                (
                    "materialProperties"
                ).subOrEmptyDict(materials[mati]).found("species")
            )
            {
                isMulticomponent = true;
            }
        }
    }
    else
    {
        isMulticomponent = (materials.size() > 1);
    }
    if (isMulticomponent && defaults.isDict("multiComponentSteady"))
    {
        if (states.found("steady"))
        {
            const_cast<dictionary&>(defaults).merge
            (
                defaults.subDict("multiComponentSteady")
            );
        }
        else if (states.found("transient"))
        {
            const_cast<dictionary&>(defaults).merge
            (
                defaults.subDict("multiComponentTransient")
            );
        }
    }
}


Xfer<dictionary> stateFunction::readDictionaries
(
    const dictionary& user,
    const stateFunction& master
) const
{

    const Time& time
    (
        refCast<const stateFunctions::globalState>(master).runTime()
    );

    dictionary existing;

    //generate const/sys sub-dicts
    existing.add("system", dictionary());
    existing.add("constant", dictionary());
    existing.add("uniform", dictionary());

    //make a copy of the stateDict and merge input
    dictionary preState(stateDict_);
    preState.merge(input_);
    preState.merge(user);

    word lclreg(regionName_);
    if (lclreg == polyMesh::defaultRegion)
    {
        lclreg = word::null;
    }

    Info<< nl << "Reusing existing dictionaries" << endl;
    Info       << "=============================" << endl;

    word dir("system");
    //read and merge system dictionaries that match pre-state entries
    if (preState.found(dir))
    {
        const dictionary& psDir(preState.subDict(dir));

        forAllConstIter(dictionary, psDir, iter)
        {
            IOobject fileDict
            (
                iter().keyword(),
                time.caseSystem(),
                lclreg,
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (fileDict.headerOk())
            {
                Info<< "Reading " << dir << "/" << iter().keyword() << endl;

                IOdictionary iodict(fileDict);

                //merge with existing
                existing.subDict(dir).add
                (
                    iter().keyword(),
                    dictionary()
                );

                existing.subDict(dir).subDict(iter().keyword()).transfer
                (
                    iodict
                );
            }
        }
    }

    dir = word("constant");
    //read and merge constant dictionaries that match pre-state entries
    if (preState.found(dir))
    {
        const dictionary& psDir(preState.subDict(dir));

        forAllConstIter(dictionary, psDir, iter)
        {
            IOobject fileDict
            (
                iter().keyword(),
                time.caseConstant(),
                lclreg,
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (fileDict.headerOk())
            {
                Info<< "Reading " << dir << "/" << iter().keyword() << endl;

                IOdictionary iodict(fileDict);

                //merge with existing
                existing.subDict(dir).add
                (
                    iter().keyword(),
                    dictionary()
                );

                existing.subDict(dir).subDict(iter().keyword()).transfer
                (
                    iodict
                );
            }
        }
    }

    dir = word("uniform");
    //read and merge constant dictionaries that match pre-state entries
    if (preState.found(dir))
    {
        const dictionary& psDir(preState.subDict(dir));

        forAllConstIter(dictionary, psDir, iter)
        {
            IOobject fileDict
            (
                iter().keyword(),
                time.timeName(),
                "uniform",
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (fileDict.headerOk())
            {
                Info<< "Reading " << dir << "/" << iter().keyword() << endl;

                IOdictionary iodict(fileDict);

                //merge with existing
                existing.subDict(dir).add
                (
                    iter().keyword(),
                    dictionary()
                );

                existing.subDict(dir).subDict(iter().keyword()).transfer
                (
                    iodict
                );
            }
        }
    }
    //Perform sync to ensure all threads are read before master write
    Pstream::barrier();

    Info<< endl;

    return existing.xfer();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stateFunction::stateFunction
(
    const dictionary& input,
    const dictionary& defaults,
    const bool distributed,
    const bool collated,
    const stateIndex& index
)
:
    regionName_(word::null),
    meshName_(regionName_),
    input_(input),
    defaults_(const_cast<dictionary&>(defaults)),
    stateId_(stateIdentifier(input, defaults)),
    modSw_(input),
    fieldNamesPtr_(),
    stateDict_(defaults.subDict("solvers").subDict("stateDictTemplate")),
    distributed_(distributed),
    collated_(collated),
    index_(index)
{}


stateFunction::stateFunction
(
    word name,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    regionName_(name),
    meshName_(meshName == word::null ? name : meshName),
    input_(),
    defaults_(const_cast<dictionary&>(defaults)),
    stateId_(stateIdentifier(input, defaults)),
    modSw_(input, master.switches()),
    fieldNamesPtr_(),
    stateDict_(defaults.subDict("solvers").subDict("stateDictTemplate")),
    distributed_(master.distributed_),
    collated_(master.collated_),
    index_(index)
{
    //assemble input components
    if (input.found("state"))
    {
        // Insert buoyancy + multi-component
        if (input.found("materials"))
        {
            insertProgramaticSettings
            (
                defaults_.subDict("solvers").subDict(stateId_),
                input
            );
        }

        // lay down base state
        input_.merge(defaults_.subDict("solvers").subDict(stateId_));

        // add module states
        //a get module names
        wordList stateList(input.lookup("state"));
        wordList moduleIds = extractModules
        (
            defaults.subDict("solvers").subDict("modules"),
            stateList
        );

        //b merge modules into data dictionary
        insertModules
        (
            defaults,
            moduleIds,
            stateId_,
            input_
        );

        //add existing dictionaries if switch true
        if (modSw_.reuseExistingDictionaries())
        {
            //loop through const/sys dictionaries, read and merge
            dictionary existing(readDictionaries(input, master));
            input_.merge(existing);
        }

        input_.merge(input);

        //merge in any state-specific defaults
        defaults_.merge(input_.subOrEmptyDict("defaults"));
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<stateFunction> stateFunction::New
(
    word name,
    const dictionary& input,
    dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
{
    word stateType = stateIdentifier(input, defaults);

    Info<< stateType << endl;

    // check if state contains "stateType" entry
    if (stateType != stateFunctions::globalState::typeName)
    {
        const dictionary& curState
        (
            defaults.subDict("solvers").subDict(stateType)
        );

        if (curState.found("stateType"))
        {
            stateType = word(curState.lookup("stateType"));
        }
    }

    const auto ctor = ctorTableLookup("function type", dictionaryConstructorTable_(), stateType);
    return autoPtr<stateFunction>
    (
        ctor(name, input, defaults, master, index, meshName)
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const dictionary& stateFunction::constant() const
{
    return stateDict_.subDict("constant");
}

dictionary& stateFunction::constant()
{
    if (!stateDict_.found("constant"))
    {
        stateDict_.add("constant", dictionary());
    }

    return stateDict_.subDict("constant");
}

const dictionary& stateFunction::system() const
{
    return stateDict_.subDict("system");
}

dictionary& stateFunction::system()
{
    if (!stateDict_.found("system"))
    {
        stateDict_.add("system", dictionary());
    }

    return stateDict_.subDict("system");
}

const dictionary& stateFunction::uniform() const
{
    return stateDict_.subDict("uniform");
}

dictionary& stateFunction::uniform()
{
    if (!stateDict_.found("uniform"))
    {
        stateDict_.add("uniform", dictionary());
    }

    return stateDict_.subDict("uniform");
}

const dictionary& stateFunction::fieldDefinitions() const
{
    return stateDict_.subDict("fields");
}

const wordList& stateFunction::fieldNames()
{
    if (!fieldNamesPtr_.valid())
    {
        fieldNamesPtr_.reset
        (
            new wordList
            (
                stateDict_.subDict("fields").toc()
            )
        );
    }

    return fieldNamesPtr_();
}


void stateFunction::initialise()
{
    //compile field components into field definition dictionaries
    dictionary& fieldMaps(stateDict_.subDict("fieldMaps"));
    dictionary& fieldDicts(stateDict_.subDict("fields"));

    if (input_.found("fieldMaps"))
    {
        fieldMaps.merge(input_.subDict("fieldMaps"));
    }

    const dictionary& fieldDB(defaults_.subDict("fields"));

    const dictionary* uiField
    (
        input_.subDictPtr("fields")
    );
    const dictionary* uiBTD
    (
        input_.subDictPtr("boundaryTypeDefaults")
    );
    const dictionary* uiBCs
    (
        input_.subDictPtr("boundaryConditions")
    );

    forAllConstIter(IDLList<entry>, fieldMaps, iter)
    {
        word targetField = iter().keyword();
        word sourceField(iter().stream());

        //start from field database entries
        fieldDicts.add(targetField, fieldDB.subDict(sourceField));

        //current field definition
        dictionary& cFieldDef(fieldDicts.subDict(targetField));

        //merge user field specifications if present
        if (uiField && uiField->found(targetField))
        {
            cFieldDef.merge
            (
                uiField->subDict(targetField)
            );
        }

        //merge user boundaryTypeDefaults
        dictionary& dictBTD
        (
            cFieldDef.subDict("fieldDefinition").subDict("boundaryTypeDefaults")
        );

        if (uiBTD)
        {
            forAllConstIter(dictionary, *uiBTD, iter)
            {
                if (iter().dict().found(targetField))
                {
                    dictBTD.remove(iter().keyword());
                    dictBTD.add
                    (
                        iter().keyword(),
                        iter().dict().subDict(targetField)
                    );
                }
            }
        }

        //add explicit boundary definitions
        dictionary& dictBCs
        (
            cFieldDef.subDict("fieldDefinition")
            .subDict("boundaryConditions")
        );

        if (uiBCs)
        {
            forAllConstIter(dictionary, *uiBCs, iter)
            {
                if (iter().dict().found(targetField))
                {
                    dictBCs.remove(iter().keyword());
                    dictBCs.add
                    (
                        iter().keyword(),
                        iter().dict().subDict(targetField)
                    );
                }
            }
        }
    }
}

void stateFunction::correct() //call this after specialised correct
{
    // merge function objects into controlDict

    if (input_.found("functions"))
    {
        system().subDict("controlDict")
            .subDict("functions").merge(input_.subDict("functions"));
    }


    //do regions scrubbing for mono-region functions

    dictionary& functionDict
    (
        system().subDict("controlDict").subDict("functions")
    );

    forAllIter(dictionary, functionDict, iter)
    {
        if (iter().isDict())
        {
            //functions defined in regions can only have single region effect
            if (regionName().size()) //master state has null name
            {
                iter().dict().remove("region");
                iter().dict().remove("regions");

                iter().dict().add("regions", wordList(1, regionName()));
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid format for function entry: "
                << iter().keyword() << nl
                << "Entry is not a dictionary." << exit(FatalError);
        }
    }
}

void stateFunction::finalise()
{
    // merge input for system/constant
    if (input_.found("system"))
    {
        system().merge(input_.subDict("system"));
    }
    // Transport properties need to be merged with false flag
    if (input_.found("constant"))
    {
        forAllIter(dictionary, input_.subDict("constant"), iter)
        {
            if (iter().isDict())
            {
                dictionary& suDict(iter().dict());
                word dictName = iter().keyword();
                if
                (
                    dictName=="transportProperties"
                 || dictName=="thermophysicalProperties")
                {
                    if (!constant().found(dictName))
                    {
                        constant().add(dictName, suDict);
                    }
                    else
                    {
                        constant().subDict(dictName).merge(suDict, false);
                    }
                }
                else
                {
                    if (!constant().found(dictName))
                    {
                        constant().add(dictName, suDict);
                    }
                    else
                    {
                        constant().subDict(dictName).merge(suDict);
                    }
                }
            }
        }
    }
    if (input_.found("uniform"))
    {
        uniform().merge(input_.subDict("uniform"), false);
    }

    if (system().found("controlDict"))
    {
        //remove and re-add functions so they sit at the end of the dictionary
        //(cosmetic)
        dictionary fobs(system().subDict("controlDict").subDict("functions"));
        system().subDict("controlDict").remove("functions");
        system().subDict("controlDict").add("functions", fobs);

        // do auto-completion of function objects here so it is after
        // controlDict merge for global
        dictionary& functionDict
        (
            system().subDict("controlDict").subDict("functions")
        );

        forAllIter(dictionary, functionDict, iter)
        {
            if (iter().isDict())
            {
                dictionary& fodict(iter().dict());

                //inject defaults
                word foType(fodict.lookup("type"));

                dictionary foBase;
                if (defaults_.subDict("functions").found(foType))
                {
                    foBase.merge(defaults_.subDict("functions").subDict(foType));
                }

                //check required entries
                wordList requiredEntries
                (
                    foBase.lookupOrDefault<wordList>
                    (
                        "requiredEntries", wordList(0, word::null)
                    )
                );

                forAll(requiredEntries, reI)
                {
                    if (!fodict.found(requiredEntries[reI]))
                    {
                        FatalErrorInFunction

                             << "Input error: caseSetup:" << regionName()
                             << ":function:" << word(iter().keyword()) << nl
                             << " Required entry `" << requiredEntries[reI]
                             << "` not found in functionObject of type "
                             << foType << "."

                             << exit(FatalError);
                    }
                }

                //merge over base
                foBase.merge(fodict);
                foBase.remove("requiredEntries");

                //re-assign
                fodict = foBase;
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid format for function entry: "
                    << iter().keyword() << nl
                    << "Entry is not a dictionary." << exit(FatalError);
            }
        }

    }
}

void stateFunction::mergeToState(const dictionary& dict)
{
    stateDict_.merge(dict);
}

void stateFunction::updateMaster(stateFunction& global)
{
    //update global dictionaries from local state
    dictionary& localGlobal
    (
        stateDict_.subDict("global")
    );

    //merge controlDict to global and remove from regional stateDict
    localGlobal.subDict("system").subDict("controlDict").merge
    (
        stateDict_.subDict("system").subDict("controlDict")
    );
    stateDict_.subDict("system").remove("controlDict");

    //merge local global from input
    if (input_.found("global"))
    {
        localGlobal.merge(input_.subDict("global"));
    }

    //- minor hack of to avoid merging libs
    fileNameList stateLibs;
    if (localGlobal.subDict("system").subDict("controlDict").found("libs"))
    {
        stateLibs = fileNameList
        (
            localGlobal.subDict("system").subDict("controlDict").lookup("libs")
        );
    }
    fileNameList globalLibs;
    bool globalLibsExist(false);
    if (global.system().subDict("controlDict").found("libs"))
    {
        globalLibs = fileNameList
        (
            global.system().subDict("controlDict").lookup("libs")
        );
        globalLibsExist = true;
    }
    forAll(stateLibs, lI)
    {
        if (!globalLibs.found(stateLibs[lI])) globalLibs.append(stateLibs[lI]);
    }

    global.mergeToState(localGlobal);

    //- ovewrite the libs based on manual merging
    if (globalLibsExist)
    {
        global.system().subDict("controlDict").remove("libs");
    }
    global.system().subDict("controlDict").add("libs", globalLibs);
}

void stateFunction::resetBoundarySwitch(const bool &reset)
{
    Switch& boundarySwitch = const_cast<Switch&>(modSw_.resetBoundaryFieldsRef());
    boundarySwitch = reset;
}

void stateFunction::writeAllDictionaries
(
    const Time& time, const wordList& uniqueSystemDicts, bool excludeUniqueDicts
) const
{
    //no need for synchronisation, since only master writes to file
    //(unless distributed, in which case it doesnt matter)

    if (collated_ || distributed_ || (Pstream::master() || !Pstream::parRun()))
    {
        writeSystemDictionaries(time, uniqueSystemDicts, excludeUniqueDicts);
        writeConstantDictionaries(time);
        writeUniformDictionaries(time);
    }
}

void stateFunction::writeSystemDictionaries
(
    const Time& time, const wordList& uniqueDicts, bool excludeUniqueDicts
) const
{
    if (switches().resetSystemDicts())
    {
        if (!excludeUniqueDicts)
        {
            // Write the mesh-unique dicts outside of the context sub-folder
            writeDirDicts
            (
                time, meshName(),
                time.caseSystem(),
                system(),
                wordList::null(),
                uniqueDicts
            );
        }
        if (meshName() == polyMesh::defaultRegion || meshName() == regionName())
        {
            //For multi-context write the other dicts for every solution context
            writeDirDicts
            (
                time,
                regionName_,
                time.caseSystem(),
                system(),
                uniqueDicts
            );
        }
        else
        {
            writeDirDicts
            (
                time,
                meshName()/regionName_,
                time.caseSystem(),
                system(),
                uniqueDicts
            );
        }
    }
}

void stateFunction::writeConstantDictionaries(const Time& time) const
{
    if (switches().resetConstDicts())
    {
        if (meshName() == polyMesh::defaultRegion || meshName() == regionName())
        {
            writeDirDicts(time, regionName_, time.caseConstant(), constant());
        }
        else
        {
            writeDirDicts
            (
                time, meshName()/regionName_, time.caseConstant(), constant()
            );
        }
    }
}

void stateFunction::writeUniformDictionaries(const Time& time) const
{
    if (switches().resetUniformDicts())
    {
        if (meshName() == polyMesh::defaultRegion || meshName() == regionName())
        {
            writeDirDicts
            (
                time, regionName_, time.timeName()+"/uniform", uniform()
            );
        }
        else
        {
            writeDirDicts
            (
                time,
                meshName()/regionName_,
                time.timeName()+"/uniform",
                uniform()
            );
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
