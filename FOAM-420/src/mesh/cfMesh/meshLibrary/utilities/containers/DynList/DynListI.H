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
    (c) Creative Fields, Ltd.
    (c) 2018-2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::checkAllocation() const
{
    if (capacity_ < UList<T>::size())
    {
        FatalErrorInFunction
            << "size " << UList<T>::size()
            << " is out of range 0 to " << capacity_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList()
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setCapacity(0); // Update UList reference to point to shortList_

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList(const label nElem)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setSize(nElem);

    # ifdef DEBUG
    checkAllocation();
    # endif
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList
(
    const label nElem,
    const T& val
)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setSize(nElem);
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList(const label nElem, const zero)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setSize(nElem);
    UList<T>::operator=(Zero);
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList
(
    const DynList<T, SizeMin>& list
)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setSize(list.size());
    UList<T>::deepCopy(list);
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>::DynList(const UList<T>& list)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    setSize(list.size());
    UList<T>::deepCopy(list);
}


template<class T, int SizeMin>
template<class ListType>
inline Foam::DynList<T, SizeMin>::DynList(const ListType& list)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    const label len = list.size();
    setSize(len);

    for (label i=0; i<len; ++i)
    {
        this->operator[](i) = list[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::DynList<T, SizeMin>::min_size() const
{
    return SizeMin;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynList<T, SizeMin>::capacity() const
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::setCapacity
(
    const label nElem
)
{
    const label nextFree = UList<T>::size();

    // Use long-list (heap)

    // Grow heapList_
    heapList_.setSize(nElem);

    // Update list reference
    UList<T>::shallowCopy(heapList_);

    // Update sizing
    capacity_ = UList<T>::size();
    UList<T>::setAddressableSize(nextFree);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::reserve
(
    const label nElem
)
{
    // Allocate more capacity if necessary
    if (nElem > capacity_)
    {
        const label newCapacity = max
        (
            SizeMin,
            max
            (
                nElem,
                // label(SizeInc + capacity_ * SizeMult / SizeDiv)
                label(2 * capacity_)
            )
        );

        setCapacity(newCapacity);
    }
    else
    {
        // Use long-list (heap)

        // Update list reference
        UList<T>::shallowCopy(heapList_);
    }
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::setSize(const label nElem)
{
    setCapacity(nElem);
    UList<T>::setAddressableSize(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::setSize
(
    const label nElem,
    const T& val
)
{
    label nextFree = List<T>::size();
    setSize(nElem);

    // Set new elements to constant value
    while (nextFree < nElem)
    {
        this->operator[](nextFree++) = val;
    }
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::resize(const label nElem)
{
    setSize(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::resize
(
    const label nElem,
    const T& val
)
{
    setSize(nElem, val);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::clear()
{
    UList<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::clearStorage()
{
    heapList_.clear();
    UList<T>::setAddressableSize(0);
    setCapacity(0);  // Also updates UList reference to point to shortList_
}


template<class T, int SizeMin>
inline Foam::label Foam::DynList<T, SizeMin>::expandStorage()
{
    const label nextFree = UList<T>::size();

    // Allow addressing into the entire list
    UList<T>::setAddressableSize(capacity_);

    return nextFree;
}


template<class T, int SizeMin>
inline Foam::DynList<T, SizeMin>&
Foam::DynList<T, SizeMin>::shrink()
{
    setCapacity(UList<T>::size());

    return *this;
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::append(const T& val)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    const label idx = UList<T>::size();

    if (idx >= capacity_)
    {
        const label newSize = 2*capacity_ + 2;
        setCapacity(newSize);
    }

    UList<T>::setAddressableSize(idx + 1);

    # ifdef DEBUG
    checkAllocation();
    # endif

    this->operator[](idx) = val;
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::appendUniq(const T& val)
{
    if (!UList<T>::found(val))
    {
        this->append(val);

        # ifdef DEBUG
        checkAllocation();
        # endif
    }
}


template<class T, int SizeMin>
inline T Foam::DynList<T, SizeMin>::remove()
{
    const label idx = UList<T>::size() - 1;

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = UList<T>::operator[](idx);

    UList<T>::setAddressableSize(idx);

    return val;
}


template<class T, int SizeMin>
inline T Foam::DynList<T, SizeMin>::removeElement(const label i)
{
    // same as DynamicList remove(label, true)

    UList<T>::swapLast(i);
    return this->remove();
}


template<class T, int SizeMin>
inline T& Foam::DynList<T, SizeMin>::newElmt(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator()(i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynList<T, SizeMin>::operator()(const label i)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    if (i >= UList<T>::size())
    {
        setSize(i + 1);
    }

    const label sz = Foam::max(UList<T>::size(), i+1);

    if (sz >= capacity_)
    {
        setCapacity(2 * sz+1);
    }

    # ifdef DEBUG
    checkAllocation();
    # endif

    return this->operator[](i);
}

template<class T, int staticSize>
inline const T& Foam::DynList<T, staticSize>::operator[](const label i) const
{
    # ifdef FULLDEBUG
    UList<T>::checkIndex(i);
    # endif

    return UList<T>::operator[](i);
}

template<class T, int staticSize>
inline T& Foam::DynList<T, staticSize>::operator[](const label i)
{
    # ifdef FULLDEBUG
    UList<T>::checkIndex(i);
    # endif

    return UList<T>::operator[](i);
}

template<class T, int staticSize>
inline Foam::label Foam::DynList<T, staticSize>::fcIndex
(
    const label index,
    const label offset
) const
{
    return (index == this->size()-offset ? 0 : index+1);
}

template<class T, int staticSize>
inline Foam::label Foam::DynList<T, staticSize>::rcIndex
(
    const label index,
    const label offset
) const
{
    return (index ? index-offset : this->size()-1);
}

template<class T, int staticSize>
inline const T& Foam::DynList<T, staticSize>::fcElement
(
    const label index,
    const label offset
) const
{
    return operator[](fcIndex(index, offset));
}

template<class T, int staticSize>
inline const T& Foam::DynList<T, staticSize>::rcElement
(
    const label index,
    const label offset
) const
{
    return operator[](rcIndex(index, offset));
}



template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynList<T, SizeMin>::operator=
(
    const Foam::DynList<T, SizeMin>& list
)
{
    # ifdef DEBUG
    checkAllocation();
    # endif

    const label nElem = list.size();
    setCapacity(nElem);
    UList<T>::setAddressableSize(nElem);

    # ifdef DEBUG
    checkAllocation();
    # endif

    for (label i=0; i<nElem; ++i)
    {
        this->operator[](i) = list[i];
    }
}


// template<class T, int SizeMin>
// template<class ListType>
// inline void Foam::DynList<T, SizeMin>::operator=
// (
//     const UList<T>& list
// )
// {
//     # ifdef DEBUG
//     checkAllocation();
//     # endif
//
//     const label nElem = list.size();
//     setCapacity(nElem);
//     UList<T>::setAddressableSize(nElem);
//
//     # ifdef DEBUG
//     checkAllocation();
//     # endif
//
//     for (label i=0; i<nElem; ++i)
//     {
//         this->operator[](i) = list[i];
//     }
// }


// ************************************************************************* //
