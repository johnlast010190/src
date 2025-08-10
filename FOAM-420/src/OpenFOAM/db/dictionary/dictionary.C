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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "db/error/error.H"
#include "global/JobInfo/JobInfo.H"
#include "db/dictionary/primitiveEntry/primitiveEntry.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "regExp.H"
#include "db/IOstreams/hashes/OSHA1stream.H"
#include "containers/Lists/DynamicList/DynamicList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(dictionary, 0);
    const dictionary dictionary::null;

    bool dictionary::writeOptionalEntries
    (
        debug::infoSwitch("writeOptionalEntries", 0)
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::entry* Foam::dictionary::lookupScopedSubEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    string::size_type dotPos = keyword.find('.');

    if (dotPos == string::npos)
    {
        // Non-scoped lookup
        return lookupEntryPtr(keyword, recursive, patternMatch);
    }
    else if (dotPos == 0)
    {
        // Starting with a '.' -> go up for every further '.' found
        ++dotPos;

        const dictionary* dictPtr = this;
        for
        (
            string::const_iterator it = keyword.begin()+1;
            it != keyword.end() && *it == '.';
            ++dotPos, ++it
        )
        {
            // Go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction
                (
                    *this
                )   << "No parent of current dictionary when searching for "
                    << keyword.substr(1)
                    << exit(FatalIOError);

                return nullptr;
            }
        }

        return dictPtr->lookupScopedSubEntryPtr
        (
            keyword.substr(dotPos),
            false,
            patternMatch
        );
    }
    else
    {
        // The first word
        const entry* entPtr = lookupScopedSubEntryPtr
        (
            keyword.substr(0, dotPos),
            false,
            patternMatch
        );

        if (!entPtr)
        {
            // Fall back to finding key with '.' so e.g. if keyword is
            // a.b.c.d it would try
            // a.b, a.b.c, a.b.c.d

            while (true)
            {
                dotPos = keyword.find('.', dotPos+1);

                entPtr = lookupEntryPtr
                (
                    keyword.substr(0, dotPos),
                    false,
                    patternMatch
                );

                if (dotPos == string::npos)
                {
                    // Parsed the whole word. Return entry or null.
                    return entPtr;
                }

                if (entPtr && entPtr->isDict())
                {
                    return entPtr->dict().lookupScopedSubEntryPtr
                    (
                        keyword.substr(dotPos),
                        false,
                        patternMatch
                    );
                }
            }
        }

        if (entPtr->isDict())
        {
            return entPtr->dict().lookupScopedSubEntryPtr
            (
                keyword.substr(dotPos),
                false,
                patternMatch
            );
        }
    }

    return nullptr;
}


bool Foam::dictionary::findInPatterns
(
    const bool patternMatch,
    const word& keyword,
    DLList<entry*>::const_iterator& wcLink,
    DLList<autoPtr<regExp>>::const_iterator& reLink
) const
{
    if (patternEntries_.size())
    {
        while (wcLink != patternEntries_.end())
        {
            if
            (
                patternMatch
              ? reLink()->match(keyword)
              : wcLink()->keyword() == keyword
            )
            {
                return true;
            }

            ++reLink;
            ++wcLink;
        }
    }

    return false;
}


bool Foam::dictionary::findInPatterns
(
    const bool patternMatch,
    const word& keyword,
    DLList<entry*>::iterator& wcLink,
    DLList<autoPtr<regExp>>::iterator& reLink
)
{
    if (patternEntries_.size())
    {
        while (wcLink != patternEntries_.end())
        {
            if
            (
                patternMatch
              ? reLink()->match(keyword)
              : wcLink()->keyword() == keyword
            )
            {
                return true;
            }

            ++reLink;
            ++wcLink;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary(const fileName& name)
:
    dictionaryName(name),
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    dictionaryName(dict.name()),
    IDLList<entry>(dict, *this),
    parent_(parentDict)
{
    forAllIter(IDLList<entry>, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patternEntries_.insert(&iter());
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary
(
    const dictionary& dict
)
:
    dictionaryName(dict.name()),
    IDLList<entry>(dict, *this),
    parent_(dictionary::null)
{
    forAllIter(IDLList<entry>, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patternEntries_.insert(&iter());
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary
(
    const dictionary* dictPtr
)
:
    parent_(dictionary::null)
{
    if (dictPtr)
    {
        operator=(*dictPtr);
    }
}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const Xfer<dictionary>& dict
)
:
    parent_(parentDict)
{
    transfer(dict());
    name() = parentDict.name() + '.' + name();
}


Foam::dictionary::dictionary
(
    const Xfer<dictionary>& dict
)
:
    parent_(dictionary::null)
{
    transfer(dict());
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>(new dictionary(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::dictionary::topDict() const
{
    const dictionary& p = parent();

    if (&p != this && !p.name().empty())
    {
        return p.topDict();
    }
    else
    {
        return *this;
    }
}


Foam::label Foam::dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::SHA1Digest Foam::dictionary::digest() const
{
    OSHA1stream os;

    // Process entries
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        os << *iter;
    }

    return os.digest();
}


Foam::tokenList Foam::dictionary::tokens() const
{
    // Serialize dictionary into a string
    OStringStream os;
    write(os, false);
    IStringStream is(os.str());

    // Parse string as tokens
    DynamicList<token> tokens;
    token t;
    while (is.read(t))
    {
        tokens.append(t);
    }

    return tokenList(tokens.xfer());
}

void Foam::dictionary::checkITstream
(
    const ITstream& is,
    const word& keyword
) const
{
    if (is.nRemainingTokens())
    {
        const label remaining = is.nRemainingTokens();

        // Similar to SafeFatalIOError
        if (JobInfo::constructed)
        {
            OSstream& err =
                FatalIOError
                (
                    "",                 // functionName
                    "",                 // sourceFileName
                    0,                  // sourceFileLineNumber
                    this->name(),       // ioFileName
                    is.lineNumber()     // ioStartLineNumber
                );

            err << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl
                << "    ";
            is.writeList(err, 0);

            err << exit(FatalIOError);
        }
        else
        {
            std::cerr
                << nl
                << "--> FOAM FATAL IO ERROR:" << nl;

            std::cerr
                << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl;

            std::cerr
                << "file: " << this->name()
                << " at line " << is.lineNumber() << '.' << nl
                << std::endl;

            ::exit(1);
        }
    }
    else if (!is.size())
    {
        // Similar to SafeFatalIOError
        if (JobInfo::constructed)
        {
            FatalIOError
            (
                "",                 // functionName
                "",                 // sourceFileName
                0,                  // sourceFileLineNumber
                this->name(),       // ioFileName
                is.lineNumber()     // ioStartLineNumber
            )
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl
                << exit(FatalIOError);
        }
        else
        {
            std::cerr
                << nl
                << "--> FOAM FATAL IO ERROR:" << nl
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl;

            std::cerr
                << "file: " << this->name()
                << " at line " << is.lineNumber() << '.' << nl
                << std::endl;

            ::exit(1);
        }
    }
}

bool Foam::dictionary::found
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if (hashedEntries_.found(keyword))
    {
        return true;
    }
    else
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::const_iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::const_iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return true;
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.found(keyword, recursive, patternMatch);
        }
        else
        {
            return false;
        }
    }
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    HashTable<entry*>::const_iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::const_iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::const_iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return wcLink();
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.lookupEntryPtr(keyword, recursive, patternMatch);
        }
        else
        {
            return nullptr;
        }
    }

    return iter();
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return wcLink();
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return const_cast<dictionary&>(parent_).lookupEntryPtr
            (
                keyword,
                recursive,
                patternMatch
            );
        }
        else
        {
            return nullptr;
        }
    }

    return iter();
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "Keyword \"" << keyword << "\" is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return *entryPtr;
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return lookupEntry(keyword, recursive, patternMatch).stream();
}


const Foam::entry* Foam::dictionary::lookupScopedEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if (keyword[0] == ':' || keyword[0] == '^')
    {
        // Go up to top level
        const dictionary* dictPtr = this;
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }

        return dictPtr->lookupScopedSubEntryPtr
        (
            keyword.substr(1),
            false,
            patternMatch
        );
    }
    else
    {
        return lookupScopedSubEntryPtr
        (
            keyword,
            recursive,
            patternMatch
        );
    }
}


bool Foam::dictionary::substituteScopedKeyword
(
    const word& keyword,
    bool mergeEntry
)
{
    const word varName = keyword(1, keyword.size()-1);

    // Lookup the variable name in the given dictionary
    const entry* ePtr = lookupScopedEntryPtr(varName, true, true);

    // If defined insert its entries into this dictionary
    if (ePtr != nullptr)
    {
        const dictionary& addDict = ePtr->dict();

        forAllConstIter(IDLList<entry>, addDict, iter)
        {
            add(iter(), mergeEntry);
        }

        return true;
    }

    return false;
}


bool Foam::dictionary::isDict(const word& keyword) const
{
    // Find non-recursive with patterns
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->isDict();
    }
    else
    {
        return false;
    }
}


const Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "Keyword \"" << keyword << "\" is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary Foam::dictionary::subOrEmptyDict
(
    const word& keyword,
    const bool mustRead
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        if (mustRead)
        {
            FatalIOErrorInFunction
            (
                *this
            )   << "Keyword \"" << keyword << "\" is undefined in dictionary "
                << name()
                << exit(FatalIOError);
        }
        return dictionary(*this, dictionary(name() + '.' + keyword));
    }
    else
    {
        return entryPtr->dict();
    }
}


const Foam::dictionary& Foam::dictionary::optionalSubDict
(
    const word& keyword
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->dict();
    }
    else
    {
        return *this;
    }
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keys(size());

    label nKeys = 0;
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        keys[nKeys++] = iter().keyword();
    }

    return keys;
}


Foam::wordList Foam::dictionary::sortedToc() const
{
    return hashedEntries_.sortedToc();
}


Foam::List<Foam::keyType> Foam::dictionary::keys(bool patterns) const
{
    List<keyType> keys(size());

    label nKeys = 0;
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        if (iter().keyword().isPattern() ? patterns : !patterns)
        {
            keys[nKeys++] = iter().keyword();
        }
    }
    keys.setSize(nKeys);

    return keys;
}


bool Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find
    (
        entryPtr->keyword()
    );

    if (mergeEntry && iter != hashedEntries_.end())
    {
        // Merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());
            delete entryPtr;

            return true;
        }
        else
        {
            // Replace existing dictionary with entry or vice versa
            IDLList<entry>::replace(iter(), entryPtr);
            delete iter();
            hashedEntries_.erase(iter);

            if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
            {
                entryPtr->name() = name() + '.' + entryPtr->keyword();

                if (entryPtr->keyword().isPattern())
                {
                    patternEntries_.insert(entryPtr);
                    patternRegexps_.insert
                    (
                        autoPtr<regExp>(new regExp(entryPtr->keyword()))
                    );
                }

                return true;
            }
            else
            {
                IOWarningInFunction((*this))
                    << "problem replacing entry "<< entryPtr->keyword()
                    << " in dictionary " << name() << endl;

                IDLList<entry>::remove(entryPtr);
                delete entryPtr;
                return false;
            }
        }
    }

    if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
    {
        entryPtr->name() = name() + '.' + entryPtr->keyword();
        IDLList<entry>::append(entryPtr);

        if (entryPtr->keyword().isPattern())
        {
            patternEntries_.insert(entryPtr);
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(entryPtr->keyword()))
            );
        }

        return true;
    }
    else
    {
        //IOWarningInFunction((*this))
        //    << "attempt to add entry "<< entryPtr->keyword()
        //    << " which already exists in dictionary " << name()
        //    << endl;

        delete entryPtr;
        return false;
    }
}


void Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    add(e.clone(*this).ptr(), mergeEntry);
}


void Foam::dictionary::add(const keyType& k, const word& w, bool overwrite)
{
    add(new primitiveEntry(k, token(w)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const Foam::string& s,
    bool overwrite
)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const label l, bool overwrite)
{
    add(new primitiveEntry(k, token(l)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const scalar s, bool overwrite)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const dictionary& d,
    bool mergeEntry
)
{
    add(new dictionaryEntry(k, *this, d), mergeEntry);
}


void Foam::dictionary::set(entry* entryPtr)
{
    entry* existingPtr = lookupEntryPtr(entryPtr->keyword(), false, true);

    // Clear dictionary so merge acts like overwrite
    if (existingPtr && existingPtr->isDict())
    {
        existingPtr->dict().clear();
    }
    add(entryPtr, true);
}


void Foam::dictionary::set(const entry& e)
{
    set(e.clone(*this).ptr());
}


void Foam::dictionary::set(const keyType& k, const dictionary& d)
{
    set(new dictionaryEntry(k, *this, d));
}

bool Foam::dictionary::merge
(
    const dictionary& dict, bool mergeEntries, bool recursive
)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalIOErrorInFunction(*this)
            << "attempted merge to self for dictionary " << name()
            << abort(FatalIOError);
    }

    bool changed = false;

    forAllConstIter(IDLList<entry>, dict, iter)
    {
        // check for leading "~" indicating delete directive
        // other than deleting on merge, "~" prefixed entries will behave
        // identically to a standard entry

        keyType keyword(iter().keyword());

        bool deleteEntry = keyword.string::removeStart("~");

        HashTable<entry*>::iterator fnd = hashedEntries_.find(keyword);

        if (fnd != hashedEntries_.end())
        {
            // Delete any matching entry if delete directive is true
            // regardless of type
            if (deleteEntry)
            {
                remove(keyword);
            }
            // Recursively merge sub-dictionaries
            // TODO: merge without copying
            else if (recursive && fnd()->isDict() && iter().isDict())
            {
                if (fnd()->dict().merge(iter().dict(), mergeEntries, true))
                {
                    changed = true;
                }
            }
            else
            {
                add(iter().clone(*this).ptr(), mergeEntries);
                changed = true;
            }
        }
        else if (!iter().isDict() || recursive)
        {
            // Not found - just add
            add(iter().clone(*this).ptr());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    IDLList<entry>::clear();
    hashedEntries_.clear();
    patternEntries_.clear();
    patternRegexps_.clear();
}


void Foam::dictionary::transfer(dictionary& dict)
{
    // Changing parents probably doesn't make much sense,
    // but what about the names?
    name() = dict.name();

    IDLList<entry>::transfer(dict);
    hashedEntries_.transfer(dict.hashedEntries_);
    patternEntries_.transfer(dict.patternEntries_);
    patternRegexps_.transfer(dict.patternRegexps_);
}


Foam::Xfer<Foam::dictionary> Foam::dictionary::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::ITstream& Foam::dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void Foam::dictionary::operator=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    name() = rhs.name();
    clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator+=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted addition assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        set(iter().clone(*this).ptr());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

Foam::dictionary Foam::operator+
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum += dict2;
    return sum;
}


Foam::dictionary Foam::operator|
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum |= dict2;
    return sum;
}


// ************************************************************************* //
