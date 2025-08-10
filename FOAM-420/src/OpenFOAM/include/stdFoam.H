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
    (c) 2017 OpenCFD Ltd.

Namespace
    stdFoam

Description
    Includes some global templates and macros used by OpenFOAM.

    Some of the templates are defined here correspond to useful
    std templates that are part of future C++ standards, or that
    are in a state of change. Defining them here provides some additional
    control over which definition are used within the OpenFOAM code-base.

SeeAlso
    - http://en.cppreference.com/w/cpp/iterator/end
    - http://en.cppreference.com/w/cpp/iterator/begin

\*---------------------------------------------------------------------------*/

#ifndef StdFoam_H
#define StdFoam_H

#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compile-time warning for use of deprecated methods (compiler-dependent).
// Use within the class declaration.
#if defined(WIN64) || defined(WIN32)
# define FOAM_DEPRECATED(since) __attribute__((__deprecated__("Since " #since)))
# define FOAM_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))
#elif defined(__clang__)
# define FOAM_DEPRECATED(since) __attribute__((deprecated("Since " #since)))
# define FOAM_DEPRECATED_FOR(since, replacement) __attribute__((deprecated("Since " #since "; use " #replacement)))
#elif defined(__cplusplus) && (__cplusplus >= 201402L)
# define FOAM_DEPRECATED(since) [[deprecated("Since " #since)]]
# define FOAM_DEPRECATED_FOR(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
#elif defined(__GNUC__)
# define FOAM_DEPRECATED(since) __attribute__((__deprecated__("Since " #since)))
# define FOAM_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))
#else
# define FOAM_DEPRECATED(since)
# define FOAM_DEPRECATED_FOR(since, replacement)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace stdFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return iterator to the beginning of the container \a c or array.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c or array.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c or array.
//  Definition as per std::cbegin C++17
template<class C>
constexpr auto cbegin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return iterator to the end of the container \a c or array.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c or array.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(const C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c or array.
//  Definition as per std::cend C++17
template<class C>
constexpr auto cend(const C& c) -> decltype(c.end())
{
    return c.end();
}


} // End namespace stdFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef RUNTIME_DEBUG_FLAGS
template <typename T, typename U>
bool forAllConditionFn(bool condition, const T &end, const U &currentEnd, const char *file, int line)
{
    if (end != currentEnd)
    {
        fprintf(stderr, "Loop condition end changed: %s:%i\n", file, line);
        abort();
    }
    return condition;
}
#define forAllCondition(cond, end, currentEnd) \
    forAllConditionFn(cond, end, currentEnd, __FILE__, __LINE__)
#else // RUNTIME_DEBUG_FLAGS
#define forAllCondition(cond, end, currentEnd) cond
#endif // RUNTIME_DEBUG_FLAGS


//- Iterate across all elements in the \a container object of type
//  \a Container.
//  \par Usage
//  \code
//  forAllIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIters, forAllIter, forAllConstIters
#define forAllIters(container,it)                                              \
    for                                                                        \
    (                                                                          \
        auto it = stdFoam::begin(container),                                   \
        it##_end = stdFoam::end(container);                                    \
        forAllCondition(it != it##_end, it##_end, stdFoam::end(container));    \
        ++it                                                                   \
    )


//- Iterate across all elements of the \a container object with const access.
//  \par Usage
//  \code
//  forAllConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIters, forAllIter, forAllConstIter
#define forAllConstIters(container,cit)                                        \
    for                                                                        \
    (                                                                          \
        auto cit = stdFoam::cbegin(container),                                 \
        cit##_end = stdFoam::cend(container);                                  \
        forAllCondition(cit != cit##_end, cit##_end, stdFoam::cend(container));\
        ++cit                                                                  \
    )


//- Loop across all elements in \a list
// \par Usage
// \code
// forAll(anyList, i)
// {
//      statements;
// }
// \endcode
// \sa forAllReverse
#define forAll(list, i) \
    for (Foam::label i=0, i##__end = (list).size(); \
         forAllCondition(i < i##__end, i##__end, (list).size()); \
         ++i)


//- Reverse loop across all elements in \a list
//  \par Usage
//  \code
//  forAllReverse(anyList, i)
//  {
//       statements;
//  }
//  \endcode
//  \sa forAll
#ifdef RUNTIME_DEBUG_FLAGS
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1, i##__end=(list).size(); \
         forAllCondition(i>=0, i##__end, (list).size()); \
         --i)
#else
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; --i)
#endif // RUNTIME_DEBUG_FLAGS


//- Iterate across all elements in the \a container object
//  of type \a Container.
//  \par Usage
//  \code
//  forAllIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIter
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin(),                        \
        iter##__end = (container).end();                                       \
        forAllCondition(iter != iter##__end, iter##__end, (container).end());  \
        ++iter                                                                 \
    )


//- Iterate across all elements in the \a container object
//  of type \a Container with const access.
//  \par Usage
//  \code
//  forAllConstIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIter
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).cbegin(),                 \
        iter##__end = (container).cend();                                      \
        forAllCondition(iter != iter##__end, iter##__end, (container).cend()); \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
