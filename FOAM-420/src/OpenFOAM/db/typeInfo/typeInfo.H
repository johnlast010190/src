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
    (c) 2023 Esi Ltd.

Typedef
    Foam::typeInfo

Description
    Basic run-time type information using word as the type's name.
    Used to enhance the standard RTTI to cover I/O.

    The user can get the type's type name using the type info access function
    \code
        type()
    \endcode

    The reference type cast template function:
    \code
        refCast<T>(r)
    \endcode

    wraps dynamic_cast to handle the bad_cast exception and generate a
    FatalError.

    The isA function:
    \code
        isA<T>(r)
    \endcode

    returns true if r is of type T or derived from type T.

\*---------------------------------------------------------------------------*/

#ifndef typeInfo_H
#define typeInfo_H

#include "db/error/error.H"
#include "db/typeInfo/className.H"
#include <type_traits>
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// declarations (for use in header files)

// Define macros to turn off spurious warnings inside the TypeName and
// TypeNameNoDebug macros
#ifdef __clang__
#define UnsetOverrideWarning                                                   \
    _Pragma("clang diagnostic push")                                           \
    _Pragma("clang diagnostic ignored \"-Winconsistent-missing-override\"")
#define ResetOverrideWarning                                                   \
    _Pragma("clang diagnostic pop")
#else
#define UnsetOverrideWarning
#define ResetOverrideWarning
#endif

//- Declare a ClassNameNoDebug() with extra virtual type info
#define TypeNameNoDebug(TypeNameString)                                        \
    UnsetOverrideWarning                                                       \
    ClassNameNoDebug(TypeNameString);                                          \
    virtual const word& type() const { return typeName; }                      \
    ResetOverrideWarning

//- Declare a ClassName() with extra virtual type info
#define TypeName(TypeNameString)                                               \
    UnsetOverrideWarning                                                       \
    ClassName(TypeNameString);                                                 \
    virtual const word& type() const { return typeName; }                      \
    ResetOverrideWarning


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Reference type cast template function,
//  wraps dynamic_cast to handle bad_cast exception and generate a FatalError.
template<class To, class From>
inline To& dynamicCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (std::bad_cast&)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << typeid(r).name()
            << " to type " << typeid(To).name()
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


//- Reference type cast template function.
//  As per dynamicCast, but handles type names via the virtual type() method.
template<class To, class From>
inline To& refCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (std::bad_cast&)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << r.type()
            << " to type " << To::typeName
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


//- Check the typeid
template<class TestType, class Type>
inline bool isType(const Type& t)
{
    if constexpr (std::is_same_v<TestType, Type>) {
        return true;
    } else {
        return typeid(t) == typeid(TestType);
    }
}


//- Check if a dynamic_cast to typeid is possible
template<class TestType, class Type>
inline bool isA(const Type& t)
{
    const Type* tPtr = &t;
    return dynamic_cast<const TestType*>(tPtr);
}


//- Determine whether a type is a reference
template<class Type>
constexpr bool isRef = std::is_reference<Type>::value;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
