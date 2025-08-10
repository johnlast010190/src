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
    (c) 2017-2023 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(const label size)
:
    parent_type(size)
{}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable
(
    const HashPtrTable<T, Key, Hash>& ht
)
:
    parent_type(ht.capacity())
{
    for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        const T* ptr = iter.object();
        if (ptr)
        {
            this->insert(iter.key(), new T(*ptr));
        }
        else
        {
            this->insert(iter.key(), nullptr);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::~HashPtrTable()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class T, class Key, class Hash>
T* Foam::HashPtrTable<T, Key, Hash>::remove(iterator& iter)
{
    if (iter.found())
    {
        T* ptr = iter.object();
        this->parent_type::erase(iter);
        return ptr;
    }

    return nullptr;
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(iterator& iter)
{
    if (iter.found())
    {
        T* ptr = iter.object();

        if (this->parent_type::erase(iter))
        {
            if (ptr)
            {
                delete ptr;
            }

            return true;
        }
    }

    return false;
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(const Key& key)
{
    auto iter = this->find(key);
    return this->erase(iter);
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::clear()
{
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        delete iter.object();
    }

    this->parent_type::clear();
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    T* ptr
)
{
    const T* old = this->get(key);

    const bool ok = this->parent_type::set(key, ptr);

    if (ok && old != ptr)
    {
        delete const_cast<T*>(old);
    }

    return ok;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>&& ptr
)
{
    return this->set(key, ptr.release());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::operator=
(
    const HashPtrTable<T, Key, Hash>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    this->clear();

    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        const T* ptr = iter.object();
        if (ptr)
        {
            this->insert(iter.key(), ptr->clone().ptr());
        }
        else
        {
            this->insert(iter.key(), nullptr);
        }
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "containers/HashTables/HashPtrTable/HashPtrTableIO.C"

// ************************************************************************* //
