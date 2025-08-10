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

Class
    Foam::HashPtrTable

Description
    A HashTable of pointers to objects of type \<T\>.

SourceFiles
    HashPtrTable.C
    HashPtrTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashPtrTable_H
#define HashPtrTable_H

#include "containers/HashTables/HashTable/HashTable.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class T, class Key, class Hash> class HashPtrTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashPtrTable<T, Key, Hash>& tbl);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashPtrTable<T, Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                        Class HashPtrTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashPtrTable
:
    public HashTable<T*, Key, Hash>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream& is, const INew& inewt);

        //- Read from dictionary using given dictionary constructor class
        template<class INew>
        void read(const dictionary& dict, const INew& inewt);


public:

    //- The template instance used for this table
    typedef HashPtrTable<T, Key, Hash> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<T*, Key, Hash> parent_type;

    using iterator = typename parent_type::iterator;
    using const_iterator = typename parent_type::const_iterator;


    // Constructors

        //- Construct given initial table size
        HashPtrTable(const label size = 128);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        HashPtrTable(Istream& is, const INew& inewt);

        //- Construct from Istream using default Istream constructor class
        HashPtrTable(Istream& is);

        //- Construct from dictionary using given dictionary constructor class
        template<class INew>
        HashPtrTable(const dictionary& dict, const INew& inewt);

        //- Construct from dictionary with default dictionary constructor class
        HashPtrTable(const dictionary& dict);

        //- Construct as copy
        HashPtrTable(const this_type& ht);


    //- Destructor
    ~HashPtrTable();


    // Member Functions

    // Access

        //- Return const pointer associated with given entry or a nullptr
        //- if the key does not exist in the table.
        inline const T* get(const Key& key) const
        {
            // Like lookup() with a nullptr
            const const_iterator iter(this->cfind(key));
            if (iter.found())
            {
                return iter();
            }
            return nullptr;
        }

      // Edit

        //- Remove and return the pointer specified by given iterator.
        //  Includes a safeguard against the end-iterator.
        T* remove(iterator& iter);

        //- Erase an entry specified by given iterator
        //  Includes a safeguard against the end-iterator.
        bool erase(iterator& iter);

        //- Erase an entry specified by the given key
        bool erase(const Key& key);

        //- Clear all entries from table
        void clear();

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const this_type& rhs);


    // IOstream Operators

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream& is,
            HashPtrTable<T, Key, Hash>& tbl
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream& os,
            const HashPtrTable<T, Key, Hash>& tbl
        );

    // Override HashTable methods

        //- Assign a new entry, overwrites existing
        inline bool set(const Key& key, T* ptr);

        //- Assign a new entry, overwrites existing
        inline bool set(const Key& key, autoPtr<T>&& ptr);

        //- Assign a new entry, overwrites existing
        bool set(const Key& key, autoPtr<T>& ptr)
        {
            return this->set(key, std::move(ptr));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "containers/HashTables/HashPtrTable/HashPtrTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
