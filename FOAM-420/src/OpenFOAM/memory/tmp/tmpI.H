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

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include <typeinfo>

// * * * * * * * * * * * * * Private Member Operators  * * * * * * * * * * * //

template<class T>
inline void Foam::tmp<T>::operator++()
{
    ptr_->operator++();

    FOAM_ASSERT(ptr_->count() <= 1) {
        FatalErrorInFunction
            << "Attempt to create more than 2 tmps referring to"
               " the same object of type " << typeName()
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
template<class... Args>
inline Foam::tmp<T> Foam::tmp<T>::New(Args&&... args)
{
    return tmp<T>(new T(std::forward<Args>(args)...));
}


template<class T>
template<class U, class... Args>
inline Foam::tmp<T> Foam::tmp<T>::NewFrom(Args&&... args)
{
    return tmp<T>(new U(std::forward<Args>(args)...));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::tmp<T>::tmp(T* tPtr)
:
    type_(TMP),
    ptr_(tPtr)
{
    FOAM_ASSERT(!tPtr || tPtr->unique()) {
        FatalErrorInFunction
            << "Attempted construction of a " << typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const T& tRef)
:
    type_(CONST_REF),
    ptr_(const_cast<T*>(&tRef))
{}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        FOAM_ASSERT(ptr_) {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }

        operator++();
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>&& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        t.ptr_ = 0;
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool allowTransfer)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        FOAM_ASSERT(ptr_) {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }

        if (allowTransfer)
        {
            t.ptr_ = 0;
        }
        else
        {
            operator++();
        }
    }
}


template<class T>
inline Foam::tmp<T>::~tmp()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isTmp() const
{
    return type_ == TMP;
}


template<class T>
inline bool Foam::tmp<T>::empty() const
{
    return (isTmp() && !ptr_);
}


template<class T>
inline bool Foam::tmp<T>::valid() const
{
    return (!isTmp() || (isTmp() && ptr_));
}


template<class T>
inline Foam::word Foam::tmp<T>::typeName() const
{
    return "tmp<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T& Foam::tmp<T>::ref() const
{
    FOAM_ASSERT(isTmp()) {
        FatalErrorInFunction
            << "Attempt to acquire non-const reference to const object"
            << " from a " << typeName()
            << abort(FatalError);
    }

    FOAM_ASSERT(ptr_) {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp())
    {
        FOAM_ASSERT(ptr_) {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }

        FOAM_ASSERT(ptr_->unique()) {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type " << typeName()
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = 0;

        return ptr;
    }
    else
    {
        return new T(*ptr_);
    }
}


template<class T>
inline void Foam::tmp<T>::clear() const
{
    if (isTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
            ptr_ = 0;
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

#ifdef NON_CONST_TMP
template<class T>
inline T& Foam::tmp<T>::operator()()
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }

    // Const-ness is automatically cast-away which is why this operator is
    // deprecated.  Use ref() where non-const access is required.
    return *ptr_;
}
#endif


template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    FOAM_ASSERT(!isTmp() || ptr_) {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    // Return const reference
    return *ptr_;
}


template<class T>
inline Foam::tmp<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    FOAM_ASSERT(isTmp()) {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const for a " << typeName()
            << abort(FatalError);
    }

    FOAM_ASSERT(!isTmp() || ptr_) {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::tmp<T>::operator->() const
{
    FOAM_ASSERT(!isTmp() || ptr_) {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmp<T>::operator=(T* tPtr)
{
    tryReregister(tPtr, ptr_);
    clear();

    FOAM_ASSERT(tPtr) {
        FatalErrorInFunction
            << "Attempted copy of a deallocated " << typeName()
            << abort(FatalError);
    }

    FOAM_ASSERT(tPtr->unique()) {
        FatalErrorInFunction
            << "Attempted assignment of a " << typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    type_ = TMP;
    ptr_ = tPtr;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    tryReregister(t.ptr_, ptr_);
    clear();

    FOAM_ASSERT(t.isTmp()) {
        FatalErrorInFunction
            << "Attempted assignment to a const reference to an object"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }

    type_ = TMP;

    FOAM_ASSERT(t.ptr_) {
        FatalErrorInFunction
            << "Attempted assignment to a deallocated " << typeName()
            << abort(FatalError);
    }

    ptr_ = t.ptr_;
    t.ptr_ = 0;
}


inline void Foam::tryReregister(void* newPtr, void* oldPtr)
{}

// ************************************************************************* //
