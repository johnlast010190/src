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
    (c) 2016-2019 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include <typeinfo>

 // * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

 template<class T>
 template<class... Args>
 inline Foam::autoPtr<T> Foam::autoPtr<T>::New(Args&&... args)
 {
     return autoPtr<T>(new T(std::forward<Args>(args)...));
 }


 template<class T>
 template<class U, class... Args>
 inline Foam::autoPtr<T> Foam::autoPtr<T>::NewFrom(Args&&... args)
 {
     return autoPtr<T>(new U(std::forward<Args>(args)...));
 }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::autoPtr<T>::autoPtr() noexcept
:
    ptr_(nullptr)
{}


template<class T>
inline constexpr Foam::autoPtr<T>::autoPtr(std::nullptr_t) noexcept
:
    ptr_(nullptr)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(T* p) noexcept
:
    ptr_(p)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(const autoPtr<T>& ap)
:
    ptr_(ap.ptr_)
{
    ap.ptr_ = 0;
}


template<class T>
inline Foam::autoPtr<T>::autoPtr(const autoPtr<T>& ap, const bool reuse)
{
    if (reuse)
    {
        ptr_ = ap.ptr_;
        ap.ptr_ = 0;
    }
    else if (ap.valid())
    {
        ptr_ = ap().clone().ptr();
    }
    else
    {
        ptr_ = nullptr;
    }
}


template<class T>
inline Foam::autoPtr<T>::~autoPtr()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::autoPtr<T>::empty() const
{
    return !ptr_;
}


template<class T>
inline bool Foam::autoPtr<T>::valid() const
{
    return ptr_;
}


template<class T>
inline T* Foam::autoPtr<T>::ptr()
{
    T* ptr = ptr_;
    ptr_ = 0;
    return ptr;
}

template<class T>
inline T* Foam::autoPtr<T>::release() noexcept
{
    T* p = ptr_;
    ptr_ = nullptr;
    return p;
}

template<class T>
inline void Foam::autoPtr<T>::set(T* p)
{
    if (ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " already allocated"
            << abort(FatalError);
    }

    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(T* p)
{
    if (ptr_)
    {
        delete ptr_;
    }

    ptr_ = p;
}

template<class T>
inline void Foam::autoPtr<T>::reset(autoPtr<T>&& ap) noexcept
{
    reset(ap.release());
}

template<class T>
inline void Foam::autoPtr<T>::clear()
{
    reset(0);
}


template<class T>
inline T* Foam::autoPtr<T>::rawPtr()
{
    return ptr_;
}


template<class T>
inline const T* Foam::autoPtr<T>::rawPtr() const
{
    return ptr_;
}


template<class T>
inline T& Foam::autoPtr<T>::rawRef()
{
    return *ptr_;
}


template<class T>
inline const T& Foam::autoPtr<T>::rawRef() const
{
    return *ptr_;
}

template<class T>
template<class... Args>
inline Foam::autoPtr<T> Foam::autoPtr<T>::clone(Args&&... args) const
{
    if (ptr_)
    {
        return autoPtr<T>(ptr_->clone(std::forward<Args>(args)...).ptr());
    }

    return autoPtr<T>();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::autoPtr<T>::operator*()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name() << " is unallocated"
            << abort(FatalError);
    }
    return *ptr_;
}

template<class T>
inline const T& Foam::autoPtr<T>::operator*() const
{
    return const_cast<autoPtr<T>*>(this)->operator*();
}

template<class T>
inline T& Foam::autoPtr<T>::operator()()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T& Foam::autoPtr<T>::operator()() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline Foam::autoPtr<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::autoPtr<T>::operator->()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::autoPtr<T>::operator->() const
{
    return const_cast<autoPtr<T>&>(*this).operator->();
}


template<class T>
inline Foam::autoPtr<T>::operator bool() const noexcept
{
    return ptr_;
}


template<class T>
inline void Foam::autoPtr<T>::operator=(T* p)
{
    reset(p);
}


template<class T>
inline void Foam::autoPtr<T>::operator=(const autoPtr<T>& ap)
{
    if (this != &ap)
    {
        reset(const_cast<autoPtr<T>&>(ap).ptr());
    }
}


// ************************************************************************* //
