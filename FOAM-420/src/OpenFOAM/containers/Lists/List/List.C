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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "containers/Lists/List/List.H"
#include "containers/Lists/List/ListLoopM.H"
#include "containers/Lists/FixedList/FixedList.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "containers/Lists/UIndirectList/UIndirectList.H"
#include "containers/Lists/BiIndirectList/BiIndirectList.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label s)
:
    UList<T>(nullptr, s)
{
    FOAM_ASSERT(this->size_ >= 0) {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();
}


template<class T>
Foam::List<T>::List(const label s, const T& a)
:
    UList<T>(nullptr, s)
{
    FOAM_ASSERT(this->size_ >= 0) {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();

    if (this->size_)
    {
        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = a;
        List_END_FOR_ALL
    }
}


template<class T>
Foam::List<T>::List(const label s, const zero)
:
    UList<T>(nullptr, s)
{
    FOAM_ASSERT(this->size_ >= 0) {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();

    if (this->size_)
    {
        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = Zero;
        List_END_FOR_ALL
    }
}

template<class T>
Foam::List<T>::List(const UList<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    const label len = this->size_;

    if (len)
    {
        alloc();

        #ifdef USEMEMCPY
        if constexpr (contiguous<T>())
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    if (this->size_)
    {
        alloc();

        #ifdef USEMEMCPY
        if constexpr (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
template<class T2>
Foam::List<T>::List(const List<T2>& a)
:
    UList<T>(nullptr, a.size())
{
    if (this->size_)
    {
        alloc();

        List_ACCESS(T, (*this), vp);
        List_CONST_ACCESS(T2, a, ap);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = T(List_ELEM(a, ap, i));
        List_END_FOR_ALL
    }
}


template<class T>
Foam::List<T>::List(const Xfer<List<T>>& lst)
{
    transfer(lst());
}


template<class T>
Foam::List<T>::List(List<T>& a, bool reuse)
:
    UList<T>(nullptr, a.size_)
{
    if (reuse)
    {
        // Delete anything previously stored here.
        clear();

        // Harvest the organs of the provided list.
        // If the incoming list uses the local elements, we need to move them.
        this->size_ = a.size_;
        if (a.v_ == (&a.localElements[0])) {
            this->v_ = &localElements[0];

            // Move each element.
            for (label i = 0; i < a.size(); i++) {
                localElements[i] = std::move(a.localElements[i]);
            }
        } else {
            // Steal the heap allocation.
            this->v_ = a.v_;
        }

        // Not strictly necessary...
        a.v_ = nullptr;
        a.size_ = 0;
    }
    else if (this->size_)
    {
        alloc();

        #ifdef USEMEMCPY
        if constexpr (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::List<T>::List(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a, bool reuse)
:
    UList<T>(nullptr, a.size())
{
    if (reuse)
    {
        transfer(a);
    }
    else
    {
        allocCopyList(a);
    }
}


template<class T>
Foam::List<T>::List(const UList<T>& a, const labelUList& mapAddressing)
:
    UList<T>(nullptr, mapAddressing.size())
{
    if (this->size_)
    {
        // Note: cannot use List_ELEM since third argument has to be index.

        alloc();

        forAll(*this, i)
        {
            this->operator[](i) = a[mapAddressing[i]];
        }
    }
}


template<class T>
template<class InputIterator>
Foam::List<T>::List(InputIterator begIter, InputIterator endIter)
:
    List<T>(begIter, endIter, std::distance(begIter, endIter))
{}


template<class T>
template<unsigned Size>
Foam::List<T>::List(const FixedList<T, Size>& lst)
:
    UList<T>(nullptr, Size)
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(const PtrList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


// Note: using first/last is not entirely accurate.
// But since the list size is correct, last() is actually ignored.
template<class T>
Foam::List<T>::List(const SLList<T>& lst)
:
    List<T>(lst.begin(), lst.end(), lst.size())
{}


template<class T>
Foam::List<T>::List(const UIndirectList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(const BiIndirectList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(std::initializer_list<T> lst)
:
    List<T>(lst.begin(), lst.end(), lst.size())
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::~List()
{
    if (this->v_ != &localElements[0]) {
        delete[] this->v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::setSize(const label newSize)
{
    FOAM_ASSERT(newSize >= 0) {
        FatalErrorInFunction
            << "bad size " << newSize
            << abort(FatalError);
    }

    if (newSize == this->size_) {
        return;
    }

    if (newSize <= 0) {
        clear();
        return;
    }

    T* nv;
    if (newSize <= ListLocalElts<T>) {
        nv = &localElements[0];
    } else {
        nv = new T[newSize];
    }

    // If necessary, move the existing elements.
    if (this->size_ > 0 && nv != this->v_)
    {
        label toMove = min(this->size_, newSize);

        for (int i = 0; i < toMove; i++) {
            nv[i] = std::move(this->v_[i]);
        }
    }

    clear();
    this->size_ = newSize;
    this->v_ = nv;
}


template<class T>
void Foam::List<T>::setSize(const label newSize, const T& a)
{
    const label oldSize = label(this->size_);
    this->setSize(newSize);

    for (label i = oldSize; i < newSize; i++) {
        this->v_[i] = a;
    }
}


template<class T>
void Foam::List<T>::transfer(List<T>& a)
{
    clear();

    // Harvest the organs of the provided list.
    // If the incoming list uses the local elements, we need to move them.
    this->size_ = a.size_;
    if (a.v_ == (&a.localElements[0])) {
        this->v_ = &localElements[0];

        // Move each element.
        for (label i = 0; i < a.size(); i++) {
            localElements[i] = std::move(a.localElements[i]);
        }
    } else {
        // Steal the heap allocation.
        this->v_ = a.v_;
    }

    // Not strictly necessary...
    a.v_ = nullptr;
    a.size_ = 0;
}


template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
void Foam::List<T>::transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a)
{
    const label newSize = a.size_;
    transfer(static_cast<List<T>&>(a));
    this->size_ = newSize;
    a.clearStorage();
}


template<class T>
void Foam::List<T>::transfer(SortableList<T>& a)
{
    // Shrink away the sort indices
    a.shrink();
    transfer(static_cast<List<T>&>(a));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
    reAlloc(a.size_);

    if (this->size_)
    {
        #ifdef USEMEMCPY
        if constexpr (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& a)
{
    FOAM_ASSERT(this != &a) {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(static_cast<const UList<T>&>(a));
}


template<class T>
void Foam::List<T>::operator=(const SLList<T>& lst)
{
    reAlloc(lst.size());

    if (this->size_)
    {
        label i = 0;
        for (auto iter = lst.cbegin(); iter != lst.cend(); ++iter)
        {
            this->operator[](i++) = *iter;
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const UIndirectList<T>& lst)
{
    reAlloc(lst.size());
    copyList(lst);
}


template<class T>
void Foam::List<T>::operator=(const BiIndirectList<T>& lst)
{
    reAlloc(lst.size());
    copyList(lst);
}


template<class T>
void Foam::List<T>::operator=(std::initializer_list<T> lst)
{
    reAlloc(lst.size());

    auto iter = lst.begin();
    forAll(*this, i)
    {
        this->operator[](i) = *iter;
        ++iter;
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "containers/Lists/List/ListIO.C"

// ************************************************************************* //
