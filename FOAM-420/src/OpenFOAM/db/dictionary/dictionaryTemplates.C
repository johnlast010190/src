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
    (c) 2017-2021 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "db/dictionary/primitiveEntry/primitiveEntry.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::dictionary::lookup
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

    return pTraits<T>(entryPtr->stream());
}

template<class T>
T Foam::dictionary::getCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    T val;
    readCompat<T>(keyword, compat, val, matchOpt);
    return val;
}

template<class T>
bool Foam::dictionary::readCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    T& val,
    enum keyType::option matchOpt,
    bool mandatory
) const
{
    const const_searcher finder(csearchCompat(keyword, compat, matchOpt));

    if (finder.found())
    {
        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }
    else if (mandatory)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return false;
}

template<class T>
T Foam::dictionary::lookupOrDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " returning the default value '" << deflt << "'"
                << endl;
        }

        return deflt;
    }
}

template<class T>
T Foam::dictionary::lookupOrDefault
(
    const List<word>& keywords,
    const T& deflt,
    bool recursive,
    bool patternMatch
) const
{
    T value(deflt);

    bool foundMatch = false;

    forAll(keywords, ki)
    {
        const word& keyword(keywords[ki]);
        const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

        if (entryPtr)
        {
            value = pTraits<T>(entryPtr->stream());
            foundMatch = true;
            break;
        }
    }

    if (writeOptionalEntries && !foundMatch)
    {
        IOInfoInFunction(*this)
            << "None of the optional entries '" << keywords << "' are present,"
            << " returning the default value '" << deflt << "'"
            << endl;
    }

    return value;
}


template<class T>
T Foam::dictionary::lookupOrAddDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
)
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " adding and returning the default value '" << deflt << "'"
                << endl;
        }

        add(new primitiveEntry(keyword, deflt));
        return deflt;
    }
}


template<class T>
bool Foam::dictionary::readEntry
(
    const word& keyword,
    T& val,
    enum keyType::option matchOpt,
    bool readOptional
) const
{
    if (readOptional == false)
    {
        return false;
    }

    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }

    return false;
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    T& val,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        entryPtr->stream() >> val;
        return true;
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " the default value '" << val << "' will be used."
                << endl;
        }

        return false;
    }
}


template<class T>
void Foam::dictionary::add(const keyType& k, const T& t, bool overwrite)
{
    add(new primitiveEntry(k, t), overwrite);
}


template<class T>
void Foam::dictionary::set(const keyType& k, const T& t)
{
    set(new primitiveEntry(k, t));
}


// ************************************************************************* //
