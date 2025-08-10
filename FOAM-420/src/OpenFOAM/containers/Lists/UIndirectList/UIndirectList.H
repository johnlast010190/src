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
    (c) 2017-2019 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::UIndirectList

Description
    A List with indirect addressing.

    Like IndirectList but does not store addressing.

    Note the const_cast of the completeList. This is so we can use it both
    on const and non-const lists. Alternative would be to have a const_
    variant etc.

SourceFiles
    UIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef UIndirectList_H
#define UIndirectList_H

#include "containers/Lists/List/List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T> class UIndirectList;
template<class T> Ostream& operator<<(Ostream&, const UIndirectList<T>&);

/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UIndirectList
{
    // Private data

        UList<T>& completeList_;
        const labelUList& addressing_;


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline UIndirectList
        (
            const UList<T>& completeList,
            const labelUList& addr
        );

        //- Copy construct (shallow copy values and addressing)
        inline UIndirectList(const UIndirectList<T>& list)
        :
            UIndirectList<T>(list.completeList_, list.addressing_)
        {}


    // Member Functions

        // Access

            //- Return the number of elements in the list
            inline label size() const;

            //- Return true if the list is empty (ie, size() is zero).
            inline bool empty() const;

            //- Return the first element of the list.
            inline T& first();

            //- Return first element of the list.
            inline const T& first() const;

            //- Return the last element of the list.
            inline T& last();

            //- Return the last element of the list.
            inline const T& last() const;

            //- Return the forward circular index, i.e. the next index
            //  which returns to the first at the end of the list
            inline label fcIndex(const label i) const;

            //- Return the reverse circular index, i.e. the previous index
            //  which returns to the last at the beginning of the list
            inline label rcIndex(const label i) const;

            //- Return the complete list
            inline const UList<T>& completeList() const;

            //- Return the list addressing
            inline const List<label>& addressing() const;


        // Member Operators

            //- Return the addressed elements as a List
            inline List<T> operator()() const;

            //- Return non-const access to an element
            inline T& operator[](const label i);

            //- Return const access to an element
            inline const T& operator[](const label i) const;

            //- Assignment to UList of addressed elements
            inline void operator=(const UList<T>& ae);

            //- Assignment to UIndirectList of addressed elements
            inline void operator=(const UIndirectList<T>& ae);

            //- Assignment of all entries to the given value
            inline void operator=(const T& t);


    // STL type definitions

        //- Type of values the UList contains.
        typedef T value_type;

        //- Type that can be used for storing into
        //  UList::value_type objects.
        typedef T& reference;

        //- Type that can be used for storing into
        //  constant UList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //  UList iterator objects.
        typedef label difference_type;

        //- The type that can represent the size of a UList.
        typedef label size_type;


    // Writing

        //- Write the List, with line-breaks in ASCII if the list length
        //  exceeds shortListLen. Using '0' suppresses line-breaks entirely.
        //  Binary output is currently still a bit of an annoyance.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;


    // Ostream operator

        //- Write List to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <T>
        (
            Ostream& os,
            const UIndirectList<T>& L
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/Lists/UIndirectList/UIndirectListI.H"

#ifdef NoRepository
    #include "containers/Lists/UIndirectList/UIndirectListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
