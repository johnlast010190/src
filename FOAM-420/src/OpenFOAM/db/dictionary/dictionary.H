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
    (c) 2018-2023 Esi Ltd.
    (c) 2016-2022 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::dictionary

Description
    A list of keyword definitions, which are a keyword followed by a number
    of values (eg, words and numbers) or by a sub-dictionary.
    Since the dictionary format is used extensively throughout OpenFOAM for
    input/output files, there are many examples of its use.

    Dictionary keywords are a plain word or a pattern (regular expression).
    The general order for searching is as follows:
      - exact match
      - pattern match (in reverse order)
      - optional recursion into the enclosing (parent) dictionaries

    The dictionary class is the base class for IOdictionary and also serves
    as a bootstrap dictionary for the objectRegistry data dictionaries.

Note
    Within dictionaries, entries can be referenced by using the '$' syntax
    familiar from shell programming. A '.' separator is used when referencing
    sub-dictionary entries. Leading '.' prefixes can be used to specify
    an entry from a parent directory. A leading ':' or '^' prefix specifies
    starting from the top-level entry. For example,

    \verbatim
    key1        val1;
    key2        $key1;   // use key1 value from current scope
    key3        $.key1;  // use key1 value from current scope

    subdict1
    {
        key1        val1b;
        key2        $..key1; // use key1 value from parent
        subdict2
        {
            key2        val2;
            key3        $...key1; // use key1 value from grandparent
        }
    }

    key4        $^subdict1.subdict2.key3;  // lookup with absolute scoping
    \endverbatim

    It is also possible to use the '${}' syntax for clarity.

SourceFiles
    dictionary.C
    dictionaryIO.C

SeeAlso
    - Foam::entry
    - Foam::dictionaryEntry
    - Foam::primitiveEntry

\*---------------------------------------------------------------------------*/

#ifndef dictionary_H
#define dictionary_H

#include "db/dictionary/entry/entry.H"
#include "containers/LinkedLists/user/IDLList.H"
#include "containers/LinkedLists/user/DLList.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Tstreams/ITstream.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/lists/wordList.H"
#include "db/typeInfo/className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class regExp;
class dictionary;
class SHA1Digest;

Istream& operator>>(Istream&, dictionary&);
Ostream& operator<<(Ostream&, const dictionary&);

/*---------------------------------------------------------------------------*\
                        Class dictionaryName Declaration
\*---------------------------------------------------------------------------*/

//- Holds name for a dictionary
class dictionaryName
{
    // Private data

        fileName name_;

public:

    // Searching

        //- Generic const/non-const dictionary entry %searcher.
        //  A %searcher provides a uniform means of finding and returning
        //  an entry pointer as well as the dictionary \a context in which
        //  the entry was located.
        //
        //  Note that the constructors and set methods are protected such
        //  that only friends of the class can set things. This safeguards
        //  against inconsistencies in context/entry.
        template<bool Const>
        class Searcher
        {
        public:
            friend dictionary;

            //- The const/non-const type for the context and sub-dictionaries
            typedef typename std::conditional
                <Const, const dictionary, dictionary>::type dict_type;

            //- The const/non-const type for entries
            typedef typename std::conditional
                <Const, const entry, entry>::type value_type;

            //- A pointer to a const/non-const dictionary
            typedef dict_type* dict_pointer;

            //- A reference to a const/non-const dictionary
            typedef dict_type& dict_reference;

            //- A pointer to a const/non-const entry
            typedef value_type* pointer;

            //- A reference to a const/non-const entry
            typedef value_type& reference;

        protected:

            //- The dictionary context for the entry
            dict_pointer dict_;

            //- The entry or nullptr
            pointer eptr_;


            //- Construct for the given dictionary context
            Searcher(dict_pointer dict)
            :
                dict_(dict),
                eptr_(nullptr)
            {}

            //- Assign the entry
            void set(pointer eptr)
            {
                eptr_ = eptr;
            }


        public:

            //- Construct null
            Searcher()
            :
                dict_(nullptr),
                eptr_(nullptr)
            {}


            //- True if entry was found
            inline bool found() const
            {
                return eptr_;
            }

            //- The containing dictionary context
            inline dict_reference context() const
            {
                return *dict_;
            }

            //- A pointer to the entry (nullptr if not found)
            inline pointer ptr() const
            {
                return eptr_;
            }

            //- A reference to the entry (Error if not found)
            inline reference ref() const
            {
                return *eptr_;
            }

            //- True if found entry is a dictionary.
            inline bool isDict() const
            {
                return eptr_ && eptr_->isDict();
            }

            //- Pointer to the found entry as a dictionary or nullptr otherwise.
            inline dict_pointer dictPtr() const
            {
                return eptr_ && eptr_->isDict() ? eptr_->dictPtr() : nullptr;
            }

            //- Reference the found entry as a dictionary.
            //  (Error if not found, or not a dictionary).
            inline dict_reference dict() const
            {
                return eptr_->dict();
            }

            //- Permit an explicit cast to the other (const/non-const) searcher
            inline explicit operator const Searcher<!Const>&() const
            {
                return *reinterpret_cast<const Searcher<!Const>*>(this);
            }

            //- A pointer to the entry (nullptr if not found)
            pointer operator->() const
            {
                return eptr_;
            }

            //- A reference to the entry (Error if not found)
            reference operator*() const
            {
                return *eptr_;
            }
        };


        //- Searcher with const access
        typedef Searcher<true> const_searcher;

        //- Searcher with non-const access
        typedef Searcher<false> searcher;

    // Friends

        //- Declare friendship with the entry class for IO
        friend class entry;

        //- Declare friendship with the searcher classes
        friend const_searcher;
        friend searcher;

    // Constructors

        //- Construct null
        dictionaryName()
        {}

        //- Construct as copy of the given fileName
        dictionaryName(const fileName& name)
        :
            name_(name)
        {}


    // Member functions

        //- Return the dictionary name
        const fileName& name() const
        {
            return name_;
        }

        //- Return the dictionary name
        fileName& name()
        {
            return name_;
        }

        //- Return the local dictionary name (final part of scoped name)
        const word dictName() const
        {
            const word scopedName = name_.name();
            const std::string::size_type i = scopedName.rfind('.');

            if (i == std::string::npos)
            {
                return scopedName;
            }
            else
            {
                return scopedName.substr(i+1);
            }
        }
};


/*---------------------------------------------------------------------------*\
                           Class dictionary Declaration
\*---------------------------------------------------------------------------*/

class dictionary
:
    public dictionaryName,
    public IDLList<entry>
{
    // Private data

        //- Report optional keywords and values if not present in dictionary
        static bool writeOptionalEntries;

        //- HashTable of the entries held on the DL-list for quick lookup
        HashTable<entry*> hashedEntries_;

        //- Parent dictionary
        const dictionary& parent_;

        //- Entries of matching patterns
        DLList<entry*> patternEntries_;

        //- Patterns as precompiled regular expressions
        DLList<autoPtr<regExp>> patternRegexps_;

    // Typedefs

        //- The storage container
        typedef IDLList<entry> parent_type;

        typedef DLList<entry*>::iterator pattern_iterator;
        typedef DLList<entry*>::const_iterator pattern_const_iterator;

        typedef DLList<autoPtr<regExp>>::iterator regexp_iterator;
        typedef DLList<autoPtr<regExp>>::const_iterator regexp_const_iterator;

   // Private Member Functions

        //- Search using a '.' for scoping.
        //  A leading dot means to use the parent dictionary.
        //  An intermediate dot separates a sub-dictionary or sub-entry.
        //  However, the use of dots is unfortunately ambiguous.
        //  The value "a.b.c.d" could be a first-level entry, a second-level
        //  entry (eg, "a" with "b.c.d", "a.b" with "c.d" etc),
        //  or just about any other combination.
        //  The heuristic tries successively longer top-level entries
        //  until there is a suitable match.
        //
        //  \param matchOpt the search mode
        const_searcher csearchDotScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Search using a '/' for scoping.
        //  Semantics as per normal files: an intermediate "." is the current
        //  dictionary level, an intermediate ".." is the parent dictionary.
        //  Note that since a slash is not a valid word character, there is no
        //  ambiguity between separator and content.
        //  No possibility or need for recursion.
        //
        //  \param matchOpt the search mode. Recursive is ignored.
        const_searcher csearchSlashScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Find and return an entry data stream pointer if present
        //  otherwise return nullptr. Allows scoping using '.'
        const entry* lookupScopedSubEntryPtr
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        ) const;

        //- Search patterns table for exact match or regular expression match
        bool findInPatterns
        (
            const bool patternMatch,
            const word& Keyword,
            DLList<entry*>::const_iterator& wcLink,
            DLList<autoPtr<regExp>>::const_iterator& reLink
        ) const;

        //- Search patterns table for exact match or regular expression match
        bool findInPatterns
        (
            const bool patternMatch,
            const word& Keyword,
            DLList<entry*>::iterator& wcLink,
            DLList<autoPtr<regExp>>::iterator& reLink
        );


public:

    //- Declare friendship with the entry class for IO
    friend class entry;


    // Declare name of the class and its debug switch
    ClassName("dictionary");


    //- Null dictionary
    static const dictionary null;


    // Constructors

        //- Construct top-level dictionary null
        dictionary();

        //- Construct top-level empty dictionary with given name
        explicit dictionary(const fileName& name);

        //- Construct given the entry name, parent dictionary and Istream,
        //  reading entries until lastEntry or EOF
        dictionary
        (
            const fileName& name,
            const dictionary& parentDict,
            Istream& is
        );

        //- Construct top-level dictionary from Istream,
        //  reading entries until EOF
        dictionary(Istream& is);

        //- Construct top-level dictionary from Istream,
        //  reading entries until EOF, optionally keeping the header
        dictionary(Istream& is, const bool keepHeader);

        //- Construct as copy given the parent dictionary
        dictionary(const dictionary& parentDict, const dictionary&);

        //- Construct top-level dictionary as copy
        dictionary(const dictionary& dict);

        //- Construct top-level dictionary as copy from pointer to dictionary.
        //  A null pointer is treated like an empty dictionary.
        dictionary(const dictionary* dictPtr);

        //- Construct by transferring parameter contents given parent dictionary
        dictionary(const dictionary& parentDict, const Xfer<dictionary>& dict);

        //- Construct top-level dictionary by transferring parameter contents
        dictionary(const Xfer<dictionary>& dict);

        //- Construct and return clone
        autoPtr<dictionary> clone() const;

        //- Construct top-level dictionary on freestore from Istream
        static autoPtr<dictionary> New(Istream& is);


    //- Destructor
    virtual ~dictionary();


    // Member functions

        //- Return the parent dictionary
        const dictionary& parent() const
        {
            return parent_;
        }

        //- Return the top of the tree
        const dictionary& topDict() const;

        //- Return line number of first token in dictionary
        label startLineNumber() const;

        //- Return line number of last token in dictionary
        label endLineNumber() const;

        //- Return the SHA1 digest of the dictionary contents
        SHA1Digest digest() const;

        //- Return the dictionary as a list of tokens
        tokenList tokens() const;


        // Search and lookup

            //- Search dictionary for given keyword
            //  If recursive, search parent dictionaries
            //  If patternMatch, use regular expressions
            bool found
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return an entry data stream pointer if present
            //  otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions
            const entry* lookupEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream pointer for manipulation
            //  if present otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            entry* lookupEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            );

            //- Find and return an entry data stream if present otherwise error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            const entry& lookupEntry
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            ITstream& lookup
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Temporary backward compatibility
            template<class T>
            T lookupType
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const
            {
                return lookup<T>(keyword, recursive, patternMatch);
            }

            //- Find and return a T,
            //  if not found throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookup
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return a T, or return the given default value
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefault
            (
                const word& keyword,
                const T& deflt,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return a T, or return the given default value
            //  First match in keyword list is returned
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefault
            (
                const List<word>& keywords,
                const T& deflt,
                bool recursive = false,
                bool patternMatch = true
            ) const;


            //- Find and return a T, if not found return the given
            //  default value, and add to dictionary.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrAddDefault
            (
                const word& keyword,
                const T& deflt,
                bool recursive = false,
                bool patternMatch = true
            );

            //- Find entry and assign to T val.
            //- FatalIOError if it is found and the number of tokens is incorrect,
            //- or it is mandatory and not found.
            //
            //  \param keyword the keyword to search for
            //  \param val the value to read into
            //  \param matchOpt search mode (default: non-recursive with patterns)
            //  \param readOpt the entry is required/optional (default: MUST_READ)
            //
            //  \return true if the entry was read
            template<class T>
            bool readEntry
            (
                const word& keyword,
                T& val,
                enum keyType::option matchOpt = keyType::REGEX,
                bool readOpt = true
            ) const;

            //- Find an entry if present, and assign to T
            //  Returns true if the entry was found.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            bool readIfPresent
            (
                const word& keyword,
                T& val,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return an entry data stream pointer if present
            //  otherwise return nullptr. Allows scoping using '.'.
            //  Special handling for ':' at start of keyword and '..'.
            const entry* lookupScopedEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Check if entry is a sub-dictionary
            bool isDict(const word& keyword) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            const dictionary* subDictPtr(const word& keyword) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            dictionary* subDictPtr(const word& keyword);

            //- Find and return a sub-dictionary
            const dictionary& subDict(const word& keyword) const;

            //- Find and return a sub-dictionary for manipulation
            dictionary& subDict(const word& keyword);

            //- Find and return a sub-dictionary as a copy, or
            //  return an empty dictionary if the sub-dictionary does not exist
            dictionary subOrEmptyDict
            (
                const word& keyword,
                const bool mustRead = false
            ) const;

            //- Find and return a sub-dictionary if found
            //  otherwise return this dictionary
            const dictionary& optionalSubDict(const word&) const;

            //- Return the table of contents
            wordList toc() const;

            //- Return the sorted table of contents
            wordList sortedToc() const;

            //- Return the list of available keys or patterns
            List<keyType> keys(bool patterns = false) const;


        // Editing

            //- Substitute the given keyword prepended by '$' with the
            //  corresponding sub-dictionary entries
            bool substituteKeyword(const word& keyword, bool mergeEntry=false);

            //- Substitute the given scoped keyword prepended by '$' with the
            //  corresponding sub-dictionary entries
            bool substituteScopedKeyword
            (
                const word& keyword,
                bool mergeEntry=false
            );

            //- Add a new entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            bool add(entry* entryPtr, bool mergeEntry=false);

            //- Add an entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            void add(const entry& e, bool mergeEntry=false);

            //- Add a word entry
            //  optionally overwrite an existing entry
            void add(const keyType& k, const word& w, bool overwrite=false);

            //- Add a string entry
            //  optionally overwrite an existing entry
            void add(const keyType& k, const string& s, bool overwrite=false);

            //- Add a label entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const label l, bool overwrite=false);

            //- Add a scalar entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const scalar s, bool overwrite=false);

            //- Add a dictionary entry
            //  optionally merge with an existing sub-dictionary
            void add
            (
                const keyType& k,
                const dictionary& d,
                bool mergeEntry = false
            );

            //- Add a T entry
            //  optionally overwrite an existing entry
            template<class T>
            void add(const keyType& k, const T& t, bool overwrite=false);

            //- Assign a new entry, overwrite any existing entry
            void set(entry* entryPtr);

            //- Assign a new entry, overwrite any existing entry
            void set(const entry& e);

            //- Assign a dictionary entry, overwrite any existing entry
            void set(const keyType& k, const dictionary& d);

            //- Assign a T entry, overwrite any existing entry
            template<class T>
            void set(const keyType& k, const T& t);

            //- Remove an entry specified by keyword
            bool remove(const word& Keyword);

            //- Change the keyword for an entry,
            //  optionally forcing overwrite of an existing entry
            bool changeKeyword
            (
                const keyType& oldKeyword,
                const keyType& newKeyword,
                bool forceOverwrite=false
            );

            //- Merge entries from the given dictionary.
            //  Also merge sub-dictionaries as required if recursive = true
            //  EdV - added switch to prevent entry, but not dictionary merging
            //  effectively creating reverseMerge when set to false
            bool merge
            (
                const dictionary& dict,
                bool mergeEntries = true,
                bool recursive = true
            );

            //- Clear the dictionary
            void clear();

            //- Transfer the contents of the argument and annul the argument.
            void transfer(dictionary& dict);

            //- Transfer contents to the Xfer container
            Xfer<dictionary> xfer();


        // Read

        //- Check after reading if the input token stream has unconsumed
        //- tokens remaining or if there were no tokens in the first place.
        //  Emits FatalIOError
        void checkITstream(const ITstream& is, const word& keyword) const;

            //- Read dictionary from Istream
            bool read(Istream& is);

            //- Read dictionary from Istream, optionally keeping the header
            bool read(Istream& is, const bool keepHeader);


        // Write

            //- Write sub-dictionary with the keyword as its header
            void writeEntry(const keyType& keyword, Ostream& os) const;

            //- Write dictionary entries.
            //  Optionally with extra new line between entries for
            //  "top-level" dictionaries
            void writeEntries(Ostream& os, const bool extraNewLine=false) const;

            //- Write dictionary, normally with sub-dictionary formatting
            void write(Ostream& os, const bool subDict=true) const;


    // Searching

        //- Search dictionary for given keyword
        //  Default search: non-recursive with patterns.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const_searcher csearch
        (
            const word& keyword,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Search dictionary for given keyword
        //  Default search: non-recursive with patterns.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const_searcher search
        (
            const word& keyword,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Search dictionary for given keyword
        //  Default search: non-recursive with patterns.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        searcher search
        (
            const word& keyword,
            enum keyType::option = keyType::REGEX
        );

        //- Search using scoping.
        //  There are two types of scoping available:
        //  -# dot-scoping, where a '.' is used to delineate the scope
        //  -# slash-scoping, where a '/' is used to delineate the scope
        //
        //  For dot-scoping, a leading '^' traverses to the top-level
        //  dictionary, leading dots mean use the parent dictionary and an
        //  intermediate dot separates a sub-dictionary or sub-entry.
        //  However, since the use of dots is ambiguous ("a.b.c" could be
        //  an entry itself or represent a "bc" entry from dictionary "a" etc),
        //  the heuristic backtracks and attempts successively longer
        //  top-level entries until a suitable match is found.
        //
        //  For slash-scoping, semantics similar to directory structures are
        //  used. A leading '/' traverses to the top-level dictionary,
        //  a single leading or intermediate '.' references the current
        //  dictionary level. A '..' pair references the parent dictionary.
        //  Any doubled slashes are silently ignored.
        //  Since a slash is not a valid keyword character, there is no
        //  ambiguity between separator and content.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const_searcher csearchScoped
        (
            const word& keyword,
            enum keyType::option
        ) const;

        //- Search using dot or slash scoping.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const_searcher searchScoped
        (
            const word& keyword,
            enum keyType::option
        ) const;

        //- Search using dot or slash scoping.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        searcher searchScoped
        (
            const word& keyword,
            enum keyType::option
        );

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        const dictionary* cfindScopedDict(const fileName& dictPath) const;

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        const dictionary* findScopedDict(const fileName& dictPath) const;

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        dictionary* findScopedDict(const fileName& dictPath);

        //- Locate existing or create sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path could not be created
        dictionary* makeScopedDict(const fileName& dictPath);


    // Compatibility helpers

        //- Search dictionary for given keyword and any compatibility names
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const_searcher csearchCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Search dictionary for given keyword and any compatibility names
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        bool foundCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Find and return an entry pointer if present, or return a nullptr,
        //- using any compatibility names if needed.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const entry* findCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option
        ) const;

        //- Find and return an entry if present, otherwise FatalIOError,
        //- using any compatibility names if needed.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        const entry& lookupEntryCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option
        ) const;

        //- Find and return a T
        //- using any compatibility names if needed.
        //- FatalIOError if not found, or if there are excess tokens.
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        template<class T>
        T getCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Find and return an entry data stream,
        //- using any compatibility names if needed.
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        ITstream& lookupCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value
        //- using any compatibility names if needed.
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        template<class T>
        T lookupOrDefaultCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            const T& deflt,
            enum keyType::option = keyType::REGEX
        ) const;

        //- Find entry and assign to T val
        //- using any compatibility names if needed.
        //- FatalIOError if there are excess tokens.
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param val the value to read
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        //
        //  \return true if the entry was found.
        template<class T>
        bool readCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            T& val,
            enum keyType::option = keyType::REGEX,
            bool mandatory = true
        ) const;

        //- Find an entry if present, and assign to T val
        //- using any compatibility names if needed.
        //- FatalIOError if it is found and there are excess tokens.
        //  Default search: non-recursive with patterns.
        //
        //  \param compat list of old compatibility keywords and the last
        //      FOAM version for which they were used.
        //  \param val the value to read
        //  \param recursive search parent dictionaries
        //  \param patternMatch use regular expressions
        //
        //  \return true if the entry was found.
        template<class T>
        bool readIfPresentCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            T& val,
            enum keyType::option = keyType::REGEX
        ) const;


    // Member Operators
        //- Find and return entry
        ITstream& operator[](const word& keyword) const;

        void operator=(const dictionary& rhs);

        //- Include entries from the given dictionary.
        //  Warn, but do not overwrite existing entries.
        void operator+=(const dictionary& rhs);

        //- Conditionally include entries from the given dictionary.
        //  Do not overwrite existing entries.
        void operator|=(const dictionary& rhs);

        //- Unconditionally include entries from the given dictionary.
        //  Overwrite existing entries.
        void operator<<=(const dictionary& rhs);


    // IOstream operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream& is, dictionary& dict);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream& os, const dictionary& dict);
};


// Global Operators

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Warn, but do not overwrite the entries from dict1.
dictionary operator+(const dictionary& dict1, const dictionary& dict2);

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Do not overwrite the entries from dict1.
dictionary operator|(const dictionary& dict1, const dictionary& dict2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/dictionary/dictionaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
