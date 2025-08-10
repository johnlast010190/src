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
    (c) 2016-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOobjectList/IOobjectList.H"
#include "db/Time/Time.H"
#include "include/OSspecific.H"
#include "db/IOobjects/IOList/IOList.H"
#include "primitives/predicates/predicates.H"

// * * * * * * * * * * * * * * Static Functions * * * * * * * * * * * * * //

namespace Foam
{
    // Templated implementation for lookup() - file-scope
    template<class UnaryMatchPredicate>
    static IOobjectList lookupImpl
    (
        const IOobjectList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        IOobjectList results(list.size());

        forAllConstIters(list, iter)
        {
            if (matcher(iter.key()))
            {
                if (IOobject::debug)
                {
                    InfoInFunction << "Found " << iter.key() << endl;
                }

                results.insert
                (
                    iter.key(),
                    new IOobject(*(iter.object()))
                );
            }
        }

        return results;
    }


    // Templated implementation for classes() - file-scope
    template<class UnaryMatchPredicate>
    static HashTable<wordHashSet> classesImpl
    (
        const IOobjectList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        HashTable<wordHashSet> summary(2*list.size());

        // Summary (key,val) = (class-name, object-names)
        forAllConstIters(list, iter)
        {
            if (matcher(iter.key()))
            {
                // Create entry (if needed) and insert
                summary(iter.object()->headerClassName()).insert(iter.key());
            }
        }

        return summary;
    }


    // Templated implementation for names(), sortedNames() - file-scope
    template<class UnaryMatchPredicate>
    static wordList namesImpl
    (
        const IOobjectList& list,
        const word& clsName,
        const UnaryMatchPredicate& matcher,
        const bool doSort
    )
    {
        wordList objNames(list.size());

        label count = 0;
        forAllConstIters(list, iter)
        {
            if (iter()->headerClassName() == clsName && matcher(iter.key()))
            {
                objNames[count++] = iter.key();
            }
        }

        objNames.setSize(count);

        if (doSort)
        {
            Foam::sort(objNames);
        }

        return objNames;
    }
}


void Foam::IOobjectList::syncNames(wordList& objNames)
{
    if (Pstream::parRun())
    {
        // Synchronize names
        Pstream::combineGather(objNames, ListUniqueEqOp<word>());
        Pstream::combineScatter(objNames);
    }

    // Sort for consistent order on all processors
    Foam::sort(objNames);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList(const label nIoObjects)
:
    HashPtrTable<IOobject>(nIoObjects)
{}


Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobject::readOption r,
    IOobject::writeOption w,
    bool registerObject
)
:
    HashPtrTable<IOobject>()
{
    word newInstance;
    fileNameList ObjectNames = fileHandler().readObjects
    (
        db,
        instance,
        local,
        newInstance
    );

    for (const auto& objName : ObjectNames)
    {
        IOobject* objectPtr = new IOobject
        (
            objName,
            newInstance,
            local,
            db,
            r,
            w,
            registerObject
        );

        bool ok = false;
        const bool throwingIOerr = FatalIOError.throwExceptions();

        try
        {
            // Use object with local scope and current instance (no searching)
            ok = objectPtr->typeHeaderOk<IOList<label>>(false, false);
        }
        catch (Foam::IOerror& err)
        {
            Warning
                << err << nl << endl;
        }

        FatalIOError.throwExceptions(throwingIOerr);

        if (ok)
        {
            insert(objectPtr->name(), objectPtr);
        }
        else
        {
            delete objectPtr;
        }
    }
}


Foam::IOobjectList::IOobjectList(const IOobjectList& iolist)
:
    HashPtrTable<IOobject>(iolist)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOobjectList::~IOobjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(IOobject& io)
{
    return insert(io.name(), &io);
}


bool Foam::IOobjectList::remove(IOobject& io)
{
    return erase(io.name());
}


Foam::IOobject* Foam::IOobjectList::lookup(const word& name) const
{
    const_iterator iter = find(name);

    if (iter.found())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << name << endl;
        }

        return const_cast<IOobject*>(*iter);
    }
    else
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Could not find " << name << endl;
        }

        return nullptr;
    }
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRe& matcher) const
{
    return lookupImpl(*this, matcher);
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRes& matcher) const
{
    return lookupImpl(*this, matcher);
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& clsName) const
{
    IOobjectList results(size());

    forAllConstIters(*this, iter)
    {
        if (iter()->headerClassName() == clsName)
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            results.insert
            (
                iter.key(),
                new IOobject(*(iter.object()))
            );
        }
    }

    return results;
}


Foam::HashTable<Foam::wordHashSet> Foam::IOobjectList::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRe& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRes& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), false);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matcher
) const
{
    return namesImpl(*this, clsName, matcher, false);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRes& matcher
) const
{
    return namesImpl(*this, clsName, matcher, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), true);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRe& matcher
) const
{
    return namesImpl(*this, clsName, matcher, true);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRes& matcher
) const
{
    return namesImpl(*this, clsName, matcher, true);
}


Foam::wordList Foam::IOobjectList::allNames() const
{
    wordList objNames(HashPtrTable<IOobject>::toc());

    syncNames(objNames);
    return objNames;
}


const Foam::IOobject* Foam::IOobjectList::cfindObject
(
    const word& objName
) const
{
    const_iterator iter = cfind(objName);

    if (iter.found())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << objName << endl;
        }

        return iter();
    }
    else if (IOobject::debug)
    {
        InfoInFunction << "Could not find " << objName << endl;
    }

    return nullptr;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOobjectList& list)
{
    os << nl << list.size() << nl << token::BEGIN_LIST << nl;

    forAllConstIters(list, it)
    {
        os << it.key() << token::SPACE << it.object()->headerClassName() << nl;
    }

    os << token::END_LIST;
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
