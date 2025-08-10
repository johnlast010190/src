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
    (c) 2017-2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::FixedList

Description
    A 1D vector of objects of type \<T\> with a fixed size \<Size\>.

SourceFiles
    FixedList.C
    FixedListI.H
    FixedListIO.C

\*---------------------------------------------------------------------------*/

#ifndef FixedList_H
#define FixedList_H

#include "primitives/bools/bool/bool.H"
#include "primitives/ints/label/label.H"
#include "primitives/ints/uLabel/uLabel.H"
#include "primitives/hashes/Hash/Hash.H"
#include "memory/autoPtr/autoPtr.H"
#include <type_traits>
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned Size> class FixedList;

template<class T, unsigned Size>
Istream& operator>>(Istream&, FixedList<T, Size>&);

template<class T, unsigned Size>
Ostream& operator<<(Ostream&, const FixedList<T, Size>&);

template<class T> class UList;

class SLListBase;
template<class LListBase, class T> class LList;
template<class T>
using SLList = LList<SLListBase, T>;

/*---------------------------------------------------------------------------*\
                           Class FixedList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned Size>
class FixedList
{
    static_assert
    (
        Size && Size <= INT_MAX,
        "Size must be positive (non-zero) and also fit as a signed value"
    );

    // Private data

        //- Vector of values of type T of size Size.
        T v_[Size];


protected:

    // Protected Member Functions

        //- Write the FixedList with its compound type
        void writeEntry(Ostream& os) const;


public:

    //- Hashing function class.
    //  Use Hasher directly for contiguous data. Otherwise hash incrementally.
    template<class HashT=Hash<T>>
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const FixedList<T, Size>&,
            unsigned seed = 0
        ) const;
    };


    // Static Member Functions

        //- Return a null FixedList
        inline static const FixedList<T, Size>& null();


    // Constructors

        //- Default construct
        FixedList() = default;

        //- Construct from value
        explicit inline FixedList(const T& t);

        //- Construct from C-array
        explicit inline FixedList(const T lst[Size]);

        //- Construct given begin/end iterators
        //  Uses std::distance when verifying the size.
        template<class InputIterator>
        inline FixedList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list
        inline FixedList(std::initializer_list<T> lst);

        //- Construct from UList
        explicit inline FixedList(const UList<T>& lst);

        //- Construct from SLList
        explicit inline FixedList(const SLList<T>& lst);

        //- Copy constructor
        inline FixedList(const FixedList<T, Size>& lst) = default;

        //- Construct from Istream
        FixedList(Istream& is);

        //- Clone
        inline autoPtr<FixedList<T, Size>> clone() const;


    // Member Functions

        // Access

            //- Return the forward circular index, i.e. the next index
            //  which returns to the first at the end of the list
            inline label fcIndex(const label i) const;

            //- Return the reverse circular index, i.e. the previous index
            //  which returns to the last at the beginning of the list
            inline label rcIndex(const label i) const;


            //- Return a const pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline const T* cdata() const;

            //- Return a pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline T* data();

            //- Return the first element of the list
            inline T& first();

            //- Return first element of the list
            inline const T& first() const;

            //- Return the last element of the list
            inline T& last();

            //- Return the last element of the list
            inline const T& last() const;


        // Check

            //- Check start is within valid range (0 ... size-1)
            inline void checkStart(const label start) const;

            //- Check size is within valid range (0 ... size)
            inline void checkSize(const label size) const;

            //- Check index i is within valid range (0 ... size-1)
            inline void checkIndex(const label i) const;

        // Search

            //- Find index of the first occurrence of the value.
            //  Any occurrences before the start pos are ignored.
            //  Linear search.
            //  \return -1 if not found.
            label find(const T& val, label pos = 0) const;

            //- Find index of the last occurrence of the value.
            //  Any occurrences after the end pos are ignored.
            //  Linear search.
            //  \return position in list or -1 if not found.
            label rfind(const T& val, label pos = -1) const;

            //- True if the value if found in the list.
            //  Any occurrences before the start pos are ignored.
            //  Linear search.
            inline bool found(const T& val, label pos = 0) const;


        // Edit

            //- Dummy resize function
            //  needed to make FixedList consistent with List
            inline void resize(const label s);

            //- Dummy setSize function
            //  needed to make FixedList consistent with List
            inline void setSize(const label s);

            //- Copy (not transfer) the argument contents
            //  needed to make FixedList consistent with List
            void transfer(const FixedList<T, Size>& lst);


    // Member operators

        //- Return element of FixedList
        inline T& operator[](const label i);

        //- Return element of constant FixedList
        inline const T& operator[](const label i) const;

        //- Assignment to array operator. Takes linear time
        inline void operator=(const T lst[Size]);

        //- Assignment to UList operator. Takes linear time
        inline void operator=(const UList<T>& lst);

        //- Assignment to SLList operator. Takes linear time
        inline void operator=(const SLList<T>& lst);

        //- Assignment to an initializer list. Takes linear time
        inline void operator=(std::initializer_list<T> lst);

        //- Assignment of all entries to the given value
        inline void operator=(const T& t);

        //- Copy assignment
        inline FixedList<T, Size>& operator=(const FixedList<T, Size>& list) = default;


    // STL type definitions

        //- Type of values the FixedList contains
        typedef T value_type;

        //- Type that can be used for storing into
        //  FixedList::value_type objects
        typedef T& reference;

        //- Type that can be used for storing into
        //  constant FixedList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //  FixedList iterator objects
        typedef label difference_type;

        //- The type that can represent the size of a FixedList
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing FixedList
        typedef T* iterator;

        //- Return an iterator to begin traversing the FixedList
        inline iterator begin();

        //- Return an iterator to end traversing the FixedList
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing FixedList
        typedef const T* const_iterator;

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant FixedList
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant FixedList
        inline const_iterator end() const;


    // STL reverse_iterator

        //- Reverse iterator for reverse traversal of FixedList
        typedef T* reverse_iterator;

        //- Return reverse_iterator to begin reverse traversing the FixedList
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the FixedList
        inline reverse_iterator rend();


    // STL const_reverse_iterator

        //- Reverse iterator for reverse traversal of constant FixedList
        typedef const T* const_reverse_iterator;

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the FixedList
        inline label size() const;

        //- Return size of the largest possible FixedList
        inline label max_size() const;

        //- Always false since zero-sized FixedList is compile-time disabled.
        inline bool empty() const;

        //- Swap two FixedLists of the same type in constant time
        void swap(FixedList<T, Size>& a);


    // STL member operators

        //- Equality operation on FixedLists of the same type.
        //  Returns true when the FixedLists are elementwise equal
        //  (using FixedList::value_type::operator==).  Takes linear time
        bool operator==(const FixedList<T, Size>& a) const;

        //- The opposite of the equality operation. Takes linear time
        bool operator!=(const FixedList<T, Size>& a) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator<(const FixedList<T, Size>& a) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator>(const FixedList<T, Size>& a) const;

        //- Return true if !(a > b). Takes linear time
        bool operator<=(const FixedList<T, Size>& a) const;

        //- Return true if !(a < b). Takes linear time
        bool operator>=(const FixedList<T, Size>& a) const;


    // Writing

        //- Write the List as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write the List, with line-breaks in ASCII if the list length
        //  exceeds shortListLen. Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List
        friend Istream& operator>> <T, Size>
        (
            Istream& is,
            FixedList<T, Size>& L
        );

        //- Write List to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <T, Size>
        (
            Ostream& os,
            const FixedList<T, Size>& L
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/Lists/FixedList/FixedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "containers/Lists/FixedList/FixedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
