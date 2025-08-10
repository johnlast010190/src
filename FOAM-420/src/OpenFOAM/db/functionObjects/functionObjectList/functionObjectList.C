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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "db/functionObjects/functionObjectList/functionObjectList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif
#include "global/argList/argList.H"
#include "db/functionObjects/timeControl/timeControlFunctionObject.H"
//#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "primitives/Tuple2/Tuple2.H"
#include "global/etcFiles/etcFiles.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"
/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::fileName Foam::functionObjectList::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::functionObjectList::createStateDict() const
{
    // Cannot set the state dictionary on construction since Time has not
    // been fully initialised
    stateDictPtr_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "functionObjectProperties",
                time_.timeName(),
                "uniform"/word("functionObjects"),
                time_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::functionObjectList::createOutputRegistry() const
{
    objectsRegistryPtr_.reset
    (
        new objectRegistry
        (
            IOobject
            (
                "functionObjectObjects",
                time_.timeName(),
                time_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


Foam::autoPtr<Foam::functionObject> Foam::functionObjectList::remove
(
    const word& key,
    label& oldIndex
)
{
    autoPtr<functionObject> oldptr;

    auto iter = indices_.find(key);  // Index of existing functionObject

    if (iter.found())
    {
        oldIndex = *iter;

        // Remove pointer from the old list
        // oldptr = this->release(oldIndex);
        oldptr = this->set(oldIndex, 0).ptr();
        indices_.erase(iter);
    }
    else
    {
        oldIndex = -1;
    }

    return oldptr;
}


void Foam::functionObjectList::listDir
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
    // Search specified directory for functionObject configuration files
    {
        fileNameList foFiles(fileHandler().readDir(dir));
        forAll(foFiles, f)
        {
            if (foFiles[f].ext().empty())
            {
                foMap.insert(foFiles[f]);
            }
        }
    }

    // Recurse into sub-directories
    {
        fileNameList foDirs(fileHandler().readDir(dir, fileName::DIRECTORY));
        forAll(foDirs, fd)
        {
            listDir(dir/foDirs[fd], foMap);
        }
    }
}


void Foam::functionObjectList::list()
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, ed)
    {
        listDir(etcDirs[ed], foMap);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << foMap.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectList::findDict(const word& funcName)
{
    // First check if there is a functionObject dictionary file in the
    // case system directory
    fileName dictFile = stringOps::expand("$FOAM_CASE")/"system"/funcName;

    if (isFile(dictFile))
    {
        return dictFile;
    }
    else
    {
        fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

        forAll(etcDirs, i)
        {
            dictFile = search(funcName, etcDirs[i]);
            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


bool Foam::functionObjectList::readFunctionObject
(
    const string& funcNameArgs,
    dictionary& functionsDict,
    HashSet<word>& requiredFields,
    const word& region
)
{
    // Parse the optional functionObject arguments:
    //     'Q(U)' -> funcName = Q; args = (U); field = U
    //
    // Supports named arguments:
    //     'patchAverage(patch=inlet, p)' -> funcName = patchAverage;
    //         args = (patch=inlet, p); field = p

    word funcName(funcNameArgs);

    int argLevel = 0;
    wordList args;

    List<Tuple2<word, string>> namedArgs;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = funcNameArgs.begin();
        iter != funcNameArgs.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            if (argLevel == 0)
            {
                funcName = funcNameArgs.substr(start, i - start);
                start = i+1;
            }
            ++argLevel;
        }
        else if (c == ',' || c == ')')
        {
            if (argLevel == 1)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            funcNameArgs.substr(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append
                    (
                        word::validate
                        (
                            funcNameArgs.substr(start, i - start)
                        )
                    );
                }
                start = i+1;
            }

            if (c == ')')
            {
                if (argLevel == 1)
                {
                    break;
                }
                --argLevel;
            }
        }
        else if (c == '=')
        {
            argName = word::validate
            (
                funcNameArgs.substr(start, i - start)
            );

            start = i+1;
            namedArg = true;
        }

        ++i;
    }

    // Search for the functionObject dictionary
    fileName path = functionObjectList::findDict(funcName);

    if (path == fileName::null)
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    //IFstream fileStream(path);
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    dictionary funcsDict(fileStream);
    dictionary* funcDictPtr = &funcsDict;

    if (funcsDict.found(funcName) && funcsDict.isDict(funcName))
    {
        funcDictPtr = &funcsDict.subDict(funcName);
    }

    dictionary& funcDict = *funcDictPtr;

    // Insert the 'field' and/or 'fields' entry corresponding to the optional
    // arguments or read the 'field' or 'fields' entry and add the required
    // fields to requiredFields
    if (args.size() == 1)
    {
        funcDict.set("field", args[0]);
        funcDict.set("fields", args);
        requiredFields.insert(args[0]);
    }
    else if (args.size() > 1)
    {
        funcDict.set("fields", args);
        requiredFields.insert(args);
    }
    else if (funcDict.found("field"))
    {
        requiredFields.insert(word(funcDict.lookup("field")));
    }
    else if (funcDict.found("fields"))
    {
        requiredFields.insert(wordList(funcDict.lookup("fields")));
    }

    // Insert named arguments
    forAll(namedArgs, i)
    {
        IStringStream entryStream
        (
            namedArgs[i].first() + ' ' + namedArgs[i].second() + ';'
        );
        funcDict.set(entry::New(entryStream).ptr());
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    // Merge this functionObject dictionary into functionsDict
    dictionary funcArgsDict;
    funcArgsDict.add(word::validate(funcNameArgs), funcDict);
    functionsDict.merge(funcArgsDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(runTime),
    parentDict_(runTime.controlDict()),
    stateDictPtr_(),
    objectsRegistryPtr_(),
    execution_(execution),
    updated_(false)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const dictionary& parentDict,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(runTime),
    parentDict_(parentDict),
    stateDictPtr_(),
    objectsRegistryPtr_(),
    execution_(execution),
    updated_(false)
{}


Foam::autoPtr<Foam::functionObjectList> Foam::functionObjectList::New
(
    const argList& args,
    const Time& runTime,
    dictionary& controlDict,
    HashSet<word>& requiredFields
)
{
    autoPtr<functionObjectList> functionsPtr;

    controlDict.add
    (
        dictionaryEntry("functions", controlDict, dictionary::null)
    );

    dictionary& functionsDict = controlDict.subDict("functions");

    word region = word::null;

    // Set the region name if specified
    if (args.optionFound("region"))
    {
        region = args["region"];
    }

    if
    (
        args.optionFound("dict")
     || args.optionFound("func")
     || args.optionFound("funcs")
    )
    {
        if (args.optionFound("dict"))
        {
            controlDict.merge
            (
                IOdictionary
                (
                    IOobject
                    (
                        args["dict"],
                        runTime,
                        IOobject::MUST_READ_IF_MODIFIED
                    )
                )
            );
        }

        if (args.optionFound("func"))
        {
            readFunctionObject
            (
                args["func"],
                functionsDict,
                requiredFields,
                region
            );
        }

        if (args.optionFound("funcs"))
        {
            wordList funcs(args.optionLookup("funcs")());

            forAll(funcs, i)
            {
                readFunctionObject
                (
                    funcs[i],
                    functionsDict,
                    requiredFields,
                    region
                );
            }
        }

        functionsPtr.reset(new functionObjectList(runTime, controlDict));
    }
    else
    {
        functionsPtr.reset(new functionObjectList(runTime));
    }

    functionsPtr->start();

    return functionsPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::functionObjectList::triggerIndex() const
{
    label triggeri = labelMin;
    stateDict().readIfPresent("triggerIndex", triggeri);

    return triggeri;
}


void Foam::functionObjectList::resetState()
{
    // Reset (re-read) the state dictionary
    stateDictPtr_.clear();
    createStateDict();
}


Foam::IOdictionary& Foam::functionObjectList::stateDict()
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    return stateDictPtr_();
}


const Foam::IOdictionary& Foam::functionObjectList::stateDict() const
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    return stateDictPtr_();
}


Foam::objectRegistry& Foam::functionObjectList::storedObjects()
{
    if (!objectsRegistryPtr_.valid())
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
}


const Foam::objectRegistry& Foam::functionObjectList::storedObjects() const
{
    if (!objectsRegistryPtr_.valid())
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
}


void Foam::functionObjectList::clear()
{
    PtrList<functionObject>::clear();
    digests_.clear();
    indices_.clear();
    updated_ = false;
}


Foam::label Foam::functionObjectList::findObjectID(const word& name) const
{
    forAll(*this, objectI)
    {
        if (operator[](objectI).name() == name)
        {
            return objectI;
        }
    }

    return -1;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    // For safety, also force a read() when execution is turned back on
    updated_ = execution_ = false;
}


bool Foam::functionObjectList::status() const
{
    return execution_;
}


bool Foam::functionObjectList::start()
{
    return read();
}


bool Foam::functionObjectList::execute(const bool reset)
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fo, "functionObject::" + objName + "::execute");
                #endif

                ok = operator[](objectI).execute() && ok;
            }

            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fo, "functionObject::" + objName + "::write");
                #endif

                ok = operator[](objectI).write() && ok;
            }

            //reset function object after execute and write
            if (reset)
            {
                this->set(objectI, nullptr);
            }
        }
    }

    // Force writing of state dictionary after function object execution
    if (time_.outputTime())
    {
        label oldPrecision = IOstream::precision_;
        IOstream::precision_ = 16;

        stateDictPtr_->writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            time_.writeCompression(),
            true
        );

        IOstream::precision_ = oldPrecision;
    }

    return ok;
}


bool Foam::functionObjectList::execute(const label subIndex)
{
    bool ok = execution_;

    if (ok)
    {
        for (functionObject& funcObj : functions())
        {
            ok = funcObj.execute(subIndex) && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::execute
(
    const UList<wordRe>& functionNames,
    const label subIndex
)
{
    bool ok = execution_;

    if (ok && functionNames.size())
    {
        for (functionObject& funcObj : functions())
        {
            if (stringOps::match(functionNames, funcObj.name()))
            {
                ok = funcObj.execute(subIndex) && ok;
            }
        }
    }

    return ok;
}


bool Foam::functionObjectList::end()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();

            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(fo, "functionObject::" + objName + "::end");
            #endif

            ok = operator[](objectI).end() && ok;
        }
    }

    clear();

    return ok;
}


bool Foam::functionObjectList::adjustTimeStep()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();

            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(fo, "functionObject::" + objName + "::adjustTimeStep");
            #endif

            ok = operator[](objectI).adjustTimeStep() && ok;
        }
    }

    return ok;
}

Foam::label Foam::functionObjectList::functionCount
(
    const dictionary& dict
) const
{
    //number of entries created from the input dictionary
    label cSize = 0;

    //mesh regions
    wordList allRegionNames = time_.names<polyMesh>();

    forAllConstIter(dictionary, dict, iter)
    {
        const word& key = iter().keyword();

        if (!iter().isDict())
        {
            //--cSize;
        }
        else
        {
            //if the regions entry is a dictionary, skip


            // priority: regions > region > no entry
            if (dict.found("regions"))
            {
                label nMatches = 0;

                if (dict.isDict("regions"))
                //for external bc which has its own "regions" format
                {
                    ++cSize;
                    ++nMatches;
                }
                else
                {
                    wordRes rnIn(wordReList(dict.lookup("regions")));

                    forAll(allRegionNames, ari)
                    {
                        if (rnIn.match(allRegionNames[ari]))
                        {
                            // Check wildcards
                            ++cSize;
                            ++nMatches;
                        }
                    }
                }

                if (nMatches == 0)
                {
                    //no match found for regions entry
                    FatalErrorInFunction
                        << "No match found for regions list: "
                        << wordReList(dict.lookup("regions"))
                        << " specified in function "
                        << key << "."
                        << exit(FatalError);
                }
            }
            else if (dict.found("region")) //legacy specification
            {
                word regn(dict.lookup("region"));
                bool found = false;

                forAll(allRegionNames, ari)
                {
                    if (allRegionNames[ari] == regn)
                    {
                        found = true;
                        ++cSize;
                        break;
                    }
                }

                if (!found)
                {
                    //region not in mesh
                    FatalErrorInFunction
                        << "Mesh region: "
                        << word(dict.lookup("region"))
                        << " specified in function "
                        << key
                        << " not found in database."
                        << exit(FatalError);
                }
            }
            else // no specification, apply to all regions
            {
                cSize += allRegionNames.size();
            }
        }
    }

    return cSize;
}

Foam::wordList Foam::functionObjectList::functionRegions
(
    const dictionary& dict
) const
{
    wordList allRegionNames;
    HashTable<const solutionInstanceRegistry*> instanceRegistries =
        time_.lookupClass<solutionInstanceRegistry>();
    if (instanceRegistries.size())
    {
        if (dict.found("instance"))
        {
            word instanceName = dict.lookup<word>("instance");
            const solutionInstanceRegistry& solReg =
                time_.lookupObject<solutionInstanceRegistry>(instanceName);
            allRegionNames = solReg.regionNames();
        }
        else
        {
            // Use the active registry
            forAllIters(instanceRegistries, iter)
            {
                if (iter()->isActive())
                {
                    allRegionNames = iter()->regionNames();
                    break;
                }
            }
        }
    }
    else
    {
        allRegionNames = time_.names<polyMesh>();
    }

    wordList funcRegNames;

    if (dict.found("regions"))
    {
        if (dict.isDict("regions"))
        //for external bc which has its own "regions" format
        {
            funcRegNames.setSize(1);
            funcRegNames[0] = word::null;
        }
        else
        {
            wordReList regionNames(dict.lookup("regions"));
            wordRes rnIn(regionNames);

            DynamicList<word> regDymList(allRegionNames.size());

            forAll(allRegionNames, ari)
            {
                if (rnIn.match(allRegionNames[ari]))
                {
                    regDymList.append(allRegionNames[ari]);
                }
            }

            regDymList.shrink();

            funcRegNames = regDymList;
        }
    }
    else if (dict.found("region")) //legacy specification
    {
        funcRegNames.setSize(1);
        funcRegNames[0] = word(dict.lookup("region"));
    }
    else // no specification, apply to all regions
    {
        funcRegNames = allRegionNames;
    }

    return funcRegNames;
}


bool Foam::functionObjectList::read()
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    bool ok = true;
    updated_ = execution_;

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr = parentDict_.lookupEntryPtr
    (
        "functions",
        false,
        false
    );

    if (entryPtr)
    {
        PtrList<functionObject> newPtrs;
        List<SHA1Digest> newDigs;
        HashTable<label> newIndices;

        label nFunc = 0;

        #if !defined( WIN32 ) && !defined( WIN64 )
        addProfiling(fo,"functionObjects::read");
        #endif

        if (!entryPtr->isDict())
        {
            FatalIOErrorInFunction(parentDict_)
                << "'functions' entry is not a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& functionsDict = entryPtr->dict();

        const_cast<Time&>(time_).libs().open
        (
            functionsDict,
            "libs",
            functionObject::dictionaryConstructorTable_()
        );

        //set size of function list (including region duplicates)
        label newSize(functionCount(functionsDict));

        newPtrs.setSize(newSize);
        newDigs.setSize(newSize);

        forAllConstIter(dictionary, functionsDict, iter)
        {
            const word& key = iter().keyword();

            if (!iter().isDict())
            {
                if (key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry " << key << " is not a dictionary" << endl;
                }

                continue;
            }

            const dictionary& dict = iter().dict();
            bool enabled = dict.lookupOrDefault("enabled", true);

            wordList funcRegNames(functionRegions(dict));

            forAll(funcRegNames, ni)
            {
                dictionary fDict(dict);
                if (fDict.found("instance"))
                {
                    word instanceName = dict.lookup<word>("instance");
                    const solutionInstanceRegistry& solReg =
                        time_.lookupObject<solutionInstanceRegistry>(instanceName);
                    if (!solReg.isActive())
                    {
                        enabled = false;
                    }
                }
                if (funcRegNames[ni] != word::null)
                {
                    fDict.remove("regions");
                }
                fDict.remove("region");
                word newKey = key;

                if
                (
                    funcRegNames[ni] != word::null
                 && funcRegNames[ni] != polyMesh::defaultRegion
                )
                {
                    //add nothing for default region
                    fDict.add("region", funcRegNames[ni]);

                    //append region name to the function name if more
                    //than one is to be instantiated
                    if (funcRegNames.size() > 1)
                    {
                        newKey = key + word("_") + funcRegNames[ni];
                    }
                }

                newDigs[nFunc] = fDict.digest();

                label oldIndex;
                Foam::autoPtr<Foam::functionObject> objPtr = remove(newKey, oldIndex);

                if (objPtr.valid())
                {
                    if (enabled)
                    {
                        // Dictionary changed for an existing functionObject
                        if (newDigs[nFunc] != digests_[oldIndex])
                        {
                            #if !defined( WIN32 ) && !defined( WIN64 )
                            addProfiling
                            (
                                fo2,
                                "functionObject::" + objPtr->name() + "::read"
                            );
                            #endif

                            enabled = objPtr->read(fDict);
                            ok = enabled && ok;
                        }
                    }

                    if (!enabled)
                    {
                        // Delete the disabled/invalid(read) functionObject
                        // delete objPtr;
                        objPtr = nullptr;
                        continue;
                    }
                }
                else if (enabled)
                {
                    autoPtr<functionObject> foPtr;

                    // Throw FatalError, FatalIOError as exceptions
                    const bool throwingError = FatalError.throwExceptions();
                    const bool throwingIOerr = FatalIOError.throwExceptions();

                    try
                    {
                        // New functionObject
                        #if !defined( WIN32 ) && !defined( WIN64 )
                        addProfiling
                        (
                            fo2,
                            "functionObject::" + newKey + "::new"
                        );
                        #endif

                        if (functionObjects::timeControl::entriesPresent(fDict))
                        {
                            foPtr.set
                            (
                                new functionObjects::timeControl
                                (newKey, time_, fDict)
                            );
                        }
                        else
                        {
                            foPtr = functionObject::New(newKey, time_, fDict);
                        }
                    }
                    catch (Foam::IOerror& ioErr)
                    {
                        Info<< ioErr << nl << endl;
                        ::exit(1);
                    }
                    catch (Foam::error& err)
                    {
                        WarningInFunction
                            << "Caught FatalError " << err << nl << endl;
                    }

                    // Restore previous exception throwing state
                    FatalError.throwExceptions(throwingError);
                    FatalIOError.throwExceptions(throwingIOerr);

                    // If one processor only has thrown an exception (so exited the
                    // constructor) invalidate the whole functionObject
                    if (returnReduce(foPtr.valid(), andOp<bool>()))
                    {
                        objPtr = foPtr.ptr();
                    }
                    else
                    {
                        ok = false;
                    }
                }


                // Insert active functionObjects into the list
                if (objPtr.valid())
                {
                    newPtrs.set(nFunc, objPtr);
                    newIndices.insert(newKey, nFunc);
                    nFunc++;
                }

            }//region loop
        }

        newPtrs.setSize(nFunc);
        newDigs.setSize(nFunc);

        // Updating the PtrList of functionObjects deletes any
        // existing unused functionObjects
        PtrList<functionObject>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
    }
    else
    {
        PtrList<functionObject>::clear();
        digests_.clear();
        indices_.clear();
    }

    return ok;
}


bool Foam::functionObjectList::filesModified() const
{
    bool ok = false;
    if (execution_)
    {
        forAll(*this, objectI)
        {
            bool changed = operator[](objectI).filesModified();
            ok = ok || changed;
        }
    }
    return ok;
}


void Foam::functionObjectList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).updateMesh(mpm);
        }
    }
}


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).movePoints(mesh);
        }
    }
}


// ************************************************************************* //
