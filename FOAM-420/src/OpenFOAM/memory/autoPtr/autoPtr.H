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

Class
    Foam::autoPtr

Description
    An auto-pointer similar to the STL auto_ptr but with automatic casting
    to a reference to the type and with pointer allocation checking on access.

SourceFiles
    autoPtrI.H

\*---------------------------------------------------------------------------*/

#ifndef autoPtr_H
#define autoPtr_H

#include <cstddef>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtr
{
    // Public data

        //- Pointer to object
        mutable T* ptr_;


public:

    typedef T Type;


    // Factory Methods

        //- Construct autoPtr of T with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_unique, but the overload for
        //      array types is not disabled.
        template<class... Args>
        inline static autoPtr<T> New(Args&&... args);

        //- Construct autoPtr from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types.
        //  In the future check for is_convertible on the pointer types
        template<class U, class... Args>
        inline static autoPtr<T> NewFrom(Args&&... args);

    // Constructors

        //- Construct with no managed object
        inline constexpr autoPtr() noexcept;

        //- Construct with no managed object (literal nullptr)
        inline constexpr autoPtr(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer
        inline explicit autoPtr(T* p) noexcept;

        //- Construct as copy by transferring pointer to this autoPtr and
        //  setting the arguments pointer to nullptr
        inline autoPtr(const autoPtr<T>&);

        //- Construct either by transferring pointer or cloning.
        //  Should only be called with type that supports cloning
        inline autoPtr(const autoPtr<T>&, const bool reuse);


    //- Destructor, delete object if pointer is not nullptr
    inline ~autoPtr();


    // Member Functions

        // Check

            //- Return true if the autoPtr is empty (ie, no pointer set)
            inline bool empty() const;

            //- Return true if the autoPtr valid (ie, the pointer is set)
            inline bool valid() const;


        // Edit

            //- Return pointer to the managed object and release ownership.
            inline T* release() noexcept;

            //- Return object pointer for reuse
            inline T* ptr();

            //- Set pointer to that given.
            //  If object pointer already set issue a FatalError
            inline void set(T*);

            //- If object pointer already set, delete object and set to given
            //  pointer
            inline void reset(T* = 0);

            //- Delete managed object and set to new given pointer
            //  \remark This is a non-standard definition, but may provide better
            //      code documentation than a simple move assign would.
            inline void reset(autoPtr<T>&& other) noexcept;

            //- Delete object (if the pointer is valid) and set pointer to
            //  nullptr
            inline void clear();


        // Access

            //- Return the pointer, without nullptr checking.
            //  Pointer remains under autoPtr management.
            inline T* rawPtr();

            //- Const access to the pointer, without nullptr checking.
            //  Pointer remains under autoPtr management.
            inline const T* rawPtr() const;

            //- Return the reference, without nullptr checking.
            inline T& rawRef();

            //- Return the const reference, without nullptr checking.
            inline const T& rawRef() const;

        // Other

            //- Construct copy by invoking clone on underlying managed object
            //  \param args list of arguments for clone
            template<class... Args>
            inline autoPtr<T> clone(Args&&... args) const;

        // Member operators

            //- Return reference to the managed object.
            //  Fatal error if no pointer is managed
            inline T& operator*();

            //- Return const reference to the object.
            //  Fatal error if no pointer is managed
            inline const T& operator*() const;

            //- Return reference to the object data
            inline T& operator()();

            //- Return const reference to the object data
            inline const T& operator()() const;

            //- Const cast to the underlying type reference
            inline operator const T&() const;

	    //- True if the managed pointer is non-null
	    explicit inline operator bool() const noexcept;

            //- Return object pointer
            inline T* operator->();

            //- Return const object pointer
            inline const T* operator->() const;

            //- Take over the object pointer from parameter
            inline void operator=(T*);

            //- Take over the object pointer from parameter
            inline void operator=(const autoPtr<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "memory/autoPtr/autoPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
