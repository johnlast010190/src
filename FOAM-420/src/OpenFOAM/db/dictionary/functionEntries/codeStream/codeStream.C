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
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/codeStream/codeStream.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"
#include "db/IOstreams/StringStreams/StringStream.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCodeContext.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        dictionaryIstream,
        codeStream
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream,
        codeStream
    );
}
}


const Foam::word Foam::functionEntries::codeStream::codeTemplateC
    = "codeStreamTemplate.C";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable& Foam::functionEntries::codeStream::libs
(
    const dictionary& dict
)
{
    const baseIOdictionary& d = static_cast<const baseIOdictionary&>
    (
        dict.topDict()
    );
    return const_cast<Time&>(d.time()).libs();
}


bool Foam::functionEntries::codeStream::doingMasterOnlyReading
(
    const dictionary& dict
)
{
    const dictionary& topDict = dict.topDict();

    if (isA<baseIOdictionary>(topDict))
    {
        const baseIOdictionary& d = static_cast<const baseIOdictionary&>
        (
            topDict
        );

        if (debug)
        {
            Pout<< "codeStream : baseIOdictionary:" << dict.name()
                << " master-only-reading:" << d.globalObject()
                << endl;
        }

        return d.globalObject();
    }
    else
    {
        if (debug)
        {
            Pout<< "codeStream : not a baseIOdictionary:" << dict.name()
                << " master-only-reading:" << regIOobject::masterOnlyReading
                << endl;
        }

        // Fall back to regIOobject::masterOnlyReading
        return regIOobject::masterOnlyReading;
    }
}


Foam::functionEntries::codeStream::streamingFunctionType
Foam::functionEntries::codeStream::getFunction
(
    const dictionary& parentDict,
    const dictionary& codeDict
)
{
    // get code, codeInclude, codeOptions
    dynamicCodeContext context(codeDict);

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    std::string sha1Str(context.sha1().str(true));
    dynamicCode dynCode("codeStream" + sha1Str, sha1Str);

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    #if defined( WIN32 ) || defined( WIN64 )
        fileName const& libPath = dynCode.libPath();
    #else
        const fileName libPath = dynCode.libPath();
    #endif

    // see if library is loaded
    void* lib = nullptr;

    const dictionary& topDict = parentDict.topDict();

    if (isA<baseIOdictionary>(topDict))
    {
        lib = libs(parentDict).findLibrary(libPath);
    }

    if (!lib)
    {
        Info<< "Using #codeStream with " << libPath << endl;
    }


    // nothing loaded
    // avoid compilation if possible by loading an existing library
    if (!lib)
    {
        if (isA<baseIOdictionary>(topDict))
        {
            // Cached access to dl libs. Guarantees clean up upon destruction
            // of Time.
            dlLibraryTable& dlLibs = libs(parentDict);
            if (dlLibs.open(libPath, false))
            {
                lib = dlLibs.findLibrary(libPath);
            }
        }
        else
        {
            // Uncached opening of libPath. Do not complain if cannot be loaded
            lib = dlOpen(libPath, false);
        }
    }


    // create library if required
    if (!lib)
    {
        const bool create =
            Pstream::master()
         || (regIOobject::fileModificationSkew <= 0);   // not NFS

        if (create)
        {
            if (!dynCode.upToDate(context))
            {
                // filter with this context
                dynCode.reset(context);

                // compile filtered C template
                dynCode.addCompileFile(codeTemplateC);

                // define Make/options
                dynCode.setOptions(context.options());
                dynCode.setLinkLibs(context.libs());

                if (!dynCode.copyOrCreateFiles(true))
                {
                    FatalIOErrorInFunction
                    (
                        parentDict
                    )   << "Failed writing files for" << nl
                        << dynCode.libRelPath() << nl
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.makeLibso())
            {
                FatalIOErrorInFunction
                (
                    parentDict
                )   << "Failed make " << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        //- Only block if we're not doing master-only reading. (flag set by
        //  regIOobject::read, baseIOdictionary constructor)
        if
        (
           !doingMasterOnlyReading(topDict)
         && regIOobject::fileModificationSkew > 0
        )
        {
            //- Since the library has only been compiled on the master the
            //  other nodes need to pick this library up through NFS
            //  We do this by just polling a few times using the
            //  fileModificationSkew.
        #if defined( WIN32 ) || defined( WIN64 )
            off64_t mySize = fileSize(libPath,false);
        #else
            off64_t mySize = fileSize(libPath);
        #endif
            off64_t masterSize = mySize;
            Pstream::scatter(masterSize);

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }


            if (mySize < masterSize)
            {
                if (debug)
                {
                    Pout<< "Local file " << libPath
                        << " not of same size (" << mySize
                        << ") as master ("
                        << masterSize << "). Waiting for "
                        << regIOobject::fileModificationSkew
                        << " seconds." << endl;
                }
                Foam::sleep(regIOobject::fileModificationSkew);

                // Recheck local size

            #if defined(WIN32) || defined(WIN64)
                mySize = Foam::fileSize(libPath,false);
            #else
                mySize = Foam::fileSize(libPath);
            #endif

                if (mySize < masterSize)
                {
                    FatalIOErrorInFunction
                    (
                        parentDict
                    )   << "Cannot read (NFS mounted) library " << nl
                        << libPath << nl
                        << "on processor " << Pstream::myProcNo()
                        << " detected size " << mySize
                        << " whereas master size is " << masterSize
                        << " bytes." << nl
                        << "If your case is not NFS mounted"
                        << " (so distributed) set fileModificationSkew"
                        << " to 0"
                        << exit(FatalIOError);
                }
            }

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " after waiting: have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }
        }

        if (isA<baseIOdictionary>(topDict))
        {
            // Cached access to dl libs. Guarantees clean up upon destruction
            // of Time.
            dlLibraryTable& dlLibs = libs(parentDict);

            if (debug)
            {
                Pout<< "Opening cached dictionary:" << libPath << endl;
            }

            if (!dlLibs.open(libPath, false))
            {
                FatalIOErrorInFunction
                (
                    parentDict
                )   << "Failed loading library " << libPath << nl
                    << "Did you add all libraries to the 'libs' entry"
                    << " in system/controlDict?"
                    << exit(FatalIOError);
            }

            lib = dlLibs.findLibrary(libPath);
        }
        else
        {
            // Uncached opening of libPath
            if (debug)
            {
                Pout<< "Opening uncached dictionary:" << libPath << endl;
            }
            lib = dlOpen(libPath, true);
        }
    }

    bool haveLib = lib;
    if (!doingMasterOnlyReading(topDict))
    {
        reduce(haveLib, andOp<bool>());
    }

    if (!haveLib)
    {
        FatalIOErrorInFunction
        (
            parentDict
        )   << "Failed loading library " << libPath
            << " on some processors."
            << exit(FatalIOError);
    }


    // Find the function handle in the library
    streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            dlSym(lib, dynCode.codeName())
        );


    if (!function)
    {
        FatalIOErrorInFunction
        (
            parentDict
        )   << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << lib << exit(FatalIOError);
    }

    return function;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    Info<< "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        parentDict
    );

    // get code dictionary
    // must reference parent for stringOps::expand to work nicely
    dictionary codeDict("#codeStream", parentDict, is);

    streamingFunctionType function = getFunction(parentDict, codeDict);

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


bool Foam::functionEntries::codeStream::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    Info<< "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        parentDict
    );

    // get code dictionary
    // must reference parent for stringOps::expand to work nicely
    dictionary codeDict("#codeStream", parentDict, is);

    streamingFunctionType function = getFunction(parentDict, codeDict);

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    parentDict.read(resultStream);

    return true;
}


// ************************************************************************* //
