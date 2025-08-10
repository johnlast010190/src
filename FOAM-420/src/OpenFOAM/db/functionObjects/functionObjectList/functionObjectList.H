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
    (c) 2015-2019 OpenCFD Ltd.

Class
    Foam::functionObjectList

Description
    List of function objects with start(), execute() and end() functions
    that is called for each object.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectList_H
#define functionObjectList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "db/functionObjects/functionObject/functionObject.H"
#include "primitives/hashes/SHA1/SHA1Digest.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class argList;
class mapPolyMesh;
class wordRe;

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    private PtrList<functionObject>
{
    // Private data

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests
        HashTable<label> indices_;

        //- Reference to Time
        const Time& time_;

        //- The parent dictionary containing a "functions" entry
        //  This entry can either be a list or a dictionary of
        //  functionObject specifications
        const dictionary& parentDict_;

        //- Function object properties - stores state information
        mutable autoPtr<IOdictionary> stateDictPtr_;

        //- Function objects output registry
        mutable autoPtr<objectRegistry> objectsRegistryPtr_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- List of functions
        const PtrList<functionObject>& functions() const { return *this; }

        //- List of functions
        PtrList<functionObject>& functions() { return *this; }

        //- Create state dictionary - attached to Time.
        void createStateDict() const;

        //- Create registry for output objects - attached to Time.
        void createOutputRegistry() const;

        //- Remove and return the function object pointer by name,
        //- and returns the old index (into digest) via the parameter.
        //  Returns nullptr (and index -1) if it didn't exist
        autoPtr<functionObject> remove(const word& key, label& oldIndex);

        //- Search the specified directory for functionObject
        //  configuration files, add to the given map and recurse
        static void listDir(const fileName& dir, HashSet<word>& foMap);

        //- Number of function objects to be instantiated
        //  Includes multiple instantiations for multi-region cases
        label functionCount(const dictionary& dict) const;


        //- List of mesh regions a function object should be instantiated for
        wordList functionRegions(const dictionary& dict) const;

        //- No copy construct
        functionObjectList(const functionObjectList&) = delete;

        //- No copy assignment
        void operator=(const functionObjectList&) = delete;


public:

    // Static data members

        //- Default relative path ("caseDicts/postProcessing") to the
        //- directory structure containing functionObject dictionary files.
        static fileName functionObjectDictPath;


    // Constructors

        //- Construct from Time and the execution setting.
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& runTime,
            const bool execution=true
        );

        //- Construct from Time, a dictionary with "functions" entry
        //  and the execution setting.
        //  \param[in]  t - the other Time instance to construct from
        //  \param[in]  parentDict - the parent dictionary containing
        //    a "functions" entry, which can either be a list or a dictionary
        //    of functionObject specifications.
        //  \param[in]  execution - whether the function objects should execute
        //    or not. Default: true.
        functionObjectList
        (
            const Time& runTime,
            const dictionary& parentDict,
            const bool execution=true
        );

        //- Construct and return a functionObjectList for an application.
        //  If the "dict" argument is specified the functionObjectList is
        //  constructed from that dictionary which is returned as
        //  controlDict otherwise the functionObjectList is constructed
        //  from the "functions" sub-dictionary of "system/controlDict"
        static autoPtr<functionObjectList> New
        (
            const argList& args,
            const Time& runTime,
            dictionary& controlDict,
            HashSet<word>& requiredFields
        );


    //- Destructor
    ~functionObjectList() = default;


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Return the current trigger index (read from the stateDict)
        label triggerIndex() const;

        //- Reset/read state dictionary for current time
        void resetState();

        //- Write access to the state dictionary ("functionObjectProperties")
        //- registered on Time
        IOdictionary& stateDict();

        //- Const access to the state dictionary ("functionObjectProperties")
        //- registered on Time
        const IOdictionary& stateDict() const;

        //- Write access to the output objects ("functionObjectObjects")
        //- registered on Time
        objectRegistry& storedObjects();

        //- Const access to the output objects ("functionObjectObjects")
        //- registered on Time
        const objectRegistry& storedObjects() const;

        //- Clear the list of function objects
        void clear();

        //- Find the ID of a given function object by name, -1 if not found.
        label findObjectID(const word& name) const;

        //- Print a list of functionObject configuration files in
        //  user/group/shipped directories.
        //  The search scheme allows for version-specific and
        //  version-independent files using the following hierarchy:
        //  - \b user settings:
        //    - ~/.FOAMcore/\<VERSION\>/caseDicts/postProcessing
        //    - ~/.FOAMcore/caseDicts/postProcessing
        //  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is set):
        //    - $FOAM_CUSTOM_DATA_DIR/\<VERSION\>/caseDicts/postProcessing
        //    - $FOAM_CUSTOM_DATA_DIR/caseDicts/postProcessing
        //  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is not set):
        //    - $FOAM_PROJECT_DIR/../site/\<VERSION\>/caseDicts/postProcessing
        //    - $FOAM_PROJECT_DIR/../site/caseDicts/postProcessing
        //  - \b other (shipped) settings:
        //    - $FOAM_PROJECT_DIR/etc/caseDicts/postProcessing
        static void list();

        //- Search for functionObject dictionary file in
        //  user/group/shipped directories.
        //  The search scheme allows for version-specific and
        //  version-independent files using the following hierarchy:
        //  - \b user settings:
        //    - ~/.FOAMcore/\<VERSION\>/caseDicts/postProcessing
        //    - ~/.FOAMcore/caseDicts/postProcessing
        //  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is set):
        //    - $FOAM_CUSTOM_DATA_DIR/\<VERSION\>/caseDicts/postProcessing
        //    - $FOAM_CUSTOM_DATA_DIR/caseDicts/postProcessing
        //  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is not set):
        //    - $FOAM_PROJECT_DIR/../site/\<VERSION\>/caseDicts/postProcessing
        //    - $FOAM_PROJECT_DIR/../site/caseDicts/postProcessing
        //  - \b other (shipped) settings:
        //    - $FOAM_PROJECT_DIR/etc/caseDicts/postProcessing
        //
        //  \return The path of the functionObject dictionary file if found
        //  otherwise an empty path
        static fileName findDict(const word& funcName);

        //- Read the specified functionObject configuration dictionary parsing
        //  the optional arguments included in the name 'funcNameArgs0',
        //  inserting 'field' or 'fields' entries as required and merging the
        //  resulting functionObject dictionary into 'functionsDict'.  Any
        //  fields required to execute the functionObject are added to
        //  'requiredFields'
        static bool readFunctionObject
        (
            const string& funcNameArgs0,
            dictionary& functionsDict,
            HashSet<word>& requiredFields,
            const word& region = word::null
        );

        //- Read and set the function objects if their data have changed
        bool read();

        //- Switch the function objects on
        void on();

        //- Switch the function objects off
        void off();

        //- Return the execution status (on/off) of the function objects
        bool status() const;

        //- Called at the start of the time-loop
        bool start();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode).
        //  Added optional reset (used in low memory post-processing mode)
        bool execute(const bool reset = false);

        //- Execute function objects using the specified subIndex.
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const label subIndex);

        //- Execute a subset of function objects using the specified subIndex.
        //  \param functionNames names or regex of existing functions to
        //      execute
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const UList<wordRe>& functionNames, const label subIndex);

        //- Called when Time::run() determines that the time-loop exits
        bool end();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        bool adjustTimeStep();

        //- Did any file get changed during execution?
        bool filesModified() const;

        //- Update for changes of mesh
        void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
