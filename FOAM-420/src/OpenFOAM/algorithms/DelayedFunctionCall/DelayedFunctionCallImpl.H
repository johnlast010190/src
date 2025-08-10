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
    (c) 2020 Esi Ltd.

Class
    Foam::DelayedFunctionCall

Description
    Class that wraps a call to a member function with arbitrary arguments,
    so that the call can be stored and enacted later.

    This can be written far more cleanly by storing a lambda expression
    as a std::function, but unfortunately this fails in gcc 4.8 due to a bug.

    See DelayedFunctionCall for base class.

See Also
    DelyedFunctionCall.H

\*---------------------------------------------------------------------------*/

#ifndef DelayedFunctionCallImpl_H
#define DelayedFunctionCallImpl_H

#include <functional>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// This namespace wraps some messy template metaprogramming used unpack a stored
// tuple into arguments for a call to a member function pointer. This can be
// replaced by one call to std::apply() in C++17
namespace tupleToArgs
{

template<class Obj, typename MemFn, typename Tuple, bool Done, int Total, int... N>
struct callWithTupleArgsImpl
{
    static void callWithTupleArgs(Obj& o, MemFn f, Tuple&& t)
    {
        callWithTupleArgsImpl
        <Obj, MemFn, Tuple, Total == 1+sizeof...(N), Total, N..., sizeof...(N)>::
        callWithTupleArgs
        (
            o, f, std::forward<Tuple>(t)
        );
    }
};

template<class Obj, typename MemFn, typename Tuple, int Total, int... N>
struct callWithTupleArgsImpl<Obj, MemFn, Tuple, true, Total, N...>
{
    static void callWithTupleArgs(Obj& o, MemFn f, Tuple&& t)
    {
        (o.*f)(std::get<N>(std::forward<Tuple>(t))...);
    }
};

template<class Obj, typename MemFn, typename Tuple>
void callWithTupleArgs(Obj& o, MemFn f, Tuple&& t)
{
    typedef typename std::decay<Tuple>::type ttype;
    callWithTupleArgsImpl
    <
        Obj,
        MemFn,
        Tuple,
        0 == std::tuple_size<ttype>::value,
        std::tuple_size<ttype>::value
    >::
    callWithTupleArgs(o, f, std::forward<Tuple>(t));
}

}


/*---------------------------------------------------------------------------*\
                  Class DelayedFunctionCallImpl Declaration
\*---------------------------------------------------------------------------*/

// Templated implementation of the DelayedFunctionCall interface

template<class ReturnType, class ObjectType, class FuncType, typename... Targs>
class DelayedFunctionCallImpl
:
    public DelayedFunctionCall<ReturnType>
{
    ObjectType& obj_;
    FuncType func_;
    const std::tuple<Targs...> args_;

public:

    // Constructor

        DelayedFunctionCallImpl
        (
            ObjectType& obj,
            FuncType f,
            Targs&&... args
        )
        :
            // Store function pointer and args as a tuple to be called later
            obj_(obj),
            func_(f),
            args_(std::forward<Targs>(args)...)
        {}


    // Member function

        //- Call the stored member function of the stored reference to object
        void call() const override
        {
            tupleToArgs::callWithTupleArgs(obj_, func_, args_);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
