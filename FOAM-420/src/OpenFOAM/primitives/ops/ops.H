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
    (c) 2018-2019 OpenCFD Ltd.
    (c) 2021 Esi Ltd

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.

    The information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

#include "primitives/Tuple2/Tuple2.H"
#include "primitives/Pair/Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define EqOp(opName, op)                                                       \
                                                                               \
template<class T1, class T2>                                                   \
class opName##Op2                                                              \
{                                                                              \
public:                                                                        \
    __attribute__((always_inline))                                             \
    void operator()(T1& x, const T2& y) const                                  \
    {                                                                          \
        op;                                                                    \
    }                                                                          \
};                                                                             \
                                                                               \
template<class T>                                                              \
class opName##Op                                                               \
{                                                                              \
public:                                                                        \
    __attribute__((always_inline))                                             \
    void operator()(T& x, const T& y) const                                    \
    {                                                                          \
        op;                                                                    \
    }                                                                          \
};

EqOp(eq, x = y)
EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(eqSqr, x = sqr(y))
EqOp(eqExp, x = exp(y))
EqOp(eqLog, x = log(y))
EqOp(eqSinh, x = sinh(y))
EqOp(eqCosh, x = cosh(y))
EqOp(eqMag, x = mag(y))
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(maxEq, x = max(x, y))
EqOp(minEq, x = min(x, y))
EqOp(minMagSqrEq, x = (magSqr(x)<=magSqr(y) ? x : y))
EqOp(maxMagSqrEq, x = (magSqr(x)>=magSqr(y) ? x : y))
EqOp(andEq, x = (x && y))
EqOp(orEq, x  = (x || y))
EqOp(notEq, x = (x != y))

EqOp(eqMinus, x = -y)

EqOp(nopEq, (void)x)

#undef EqOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define Op(opName, op)                                                         \
                                                                               \
    template<class T, class T1, class T2>                                      \
    class opName##Op3                                                          \
    {                                                                          \
    public:                                                                    \
        [[nodiscard]] __attribute__((always_inline))                           \
        T operator()(const T1& x, const T2& y) const                           \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T1, class T2>                                               \
    class opName##Op2                                                          \
    {                                                                          \
    public:                                                                    \
        [[nodiscard]]  __attribute__((always_inline))                          \
        T1 operator()(const T1& x, const T2& y) const                          \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    class opName##Op                                                           \
    {                                                                          \
    public:                                                                    \
        [[nodiscard]] __attribute__((always_inline))                           \
        T operator()(const T& x, const T& y) const                             \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<typename>                                                         \
    struct is_ ## opName ## Op : std::false_type {};                           \
                                                                               \
    template<typename T>                                                       \
    struct is_ ## opName ## Op<opName ## Op<T>>: std::true_type {};            \
                                                                               \
    template<typename T, typename Q>                                           \
    struct is_ ## opName ## Op<opName ## Op2<T, Q>>: std::true_type {};        \
                                                                               \
    template<typename T, typename Q, typename T2>                              \
    struct is_ ## opName ## Op<opName ## Op3<T, Q, T2>>: std::true_type {};


#define weightedOp(opName, op)                                                 \
                                                                               \
    template<class Type, class CombineOp>                                      \
    class opName##WeightedOp                                                   \
    {                                                                          \
        const CombineOp& cop_;                                                 \
                                                                               \
        public:                                                                \
                                                                               \
            opName##WeightedOp(const CombineOp& cop)                           \
            :                                                                  \
                cop_(cop)                                                      \
            {}                                                                 \
\
            __attribute__((always_inline))                                     \
            void operator()                                                    \
            (                                                                  \
                Type& x,                                                       \
                const label index,                                             \
                const Type& y,                                                 \
                const scalar weight                                            \
            ) const                                                            \
            {                                                                  \
                cop_(x, op);                                                   \
            }                                                                  \
    };                                                                         \


Op(sum, x + y)

Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptPow, cmptPow(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(max, max(x, y))
Op(min, min(x, y))
Op(minMagSqr, (magSqr(x)<=magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x)>=magSqr(y) ? x : y))
Op(minMod, minMod(x, y))
Op(and, x && y)
Op(or, x || y)
Op(eqEq, x == y)
Op(less, x < y)
Op(lessEq, x <= y)
Op(greater, x > y)
Op(greaterEq, x >= y)

weightedOp(multiply, (weight*y))

#undef Op
#undef weightedOp

template<typename Op>
constexpr bool IsCommutative = false;

// Specialise to say all the things that are commutative. The T here is the
// argument type, allowing these specialisations to be type-dependent (eg.
// to express that matrix multiplication isn't commutative, but other
// flavours of multiplication might be).

template<typename T>
constexpr bool IsCommutative<sumOp<T>> = true;
// Lol, plusOp and sumOp are different types?!
template<typename T>
constexpr bool IsCommutative<plusOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<minOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<minModOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<maxOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<eqEqOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<orOp<T>> = true;
template<typename T>
constexpr bool IsCommutative<andOp<T>> = true;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Operators that carry along a label and processor number

#define EqAtOp(opName, op)                                                     \
                                                                               \
template<class T>                                                              \
class opName##AtOp                                                             \
{                                                                              \
public:                                                                        \
     __attribute__((always_inline))                                            \
    void operator()                                                            \
    (                                                                          \
        Tuple2<T,Pair<label>>& x, const Tuple2<T,Pair<label>>& y               \
    ) const                                                                    \
    {                                                                          \
        op;                                                                    \
    }                                                                          \
};

EqAtOp(minEq, x = (x.first() <= y.first() ? x : y))
EqAtOp(maxEq, x = (x.first() >= y.first() ? x : y))
EqAtOp(minMagSqrEq, x = (magSqr(x.first())<=magSqr(y.first()) ? x : y))
EqAtOp(maxMagSqrEq, x = (magSqr(x.first())>=magSqr(y.first()) ? x : y))

#undef EqAtOp


#define AtOp(opName, op)                                                       \
                                                                               \
    template<class T>                                                          \
    class opName##AtOp                                                         \
    {                                                                          \
    public:                                                                    \
        [[nodiscard]] __attribute__((always_inline))                                         \
        Tuple2<T,Pair<label>>                                                  \
        operator()                                                             \
        (                                                                      \
            const Tuple2<T,Pair<label>>& x, const Tuple2<T,Pair<label>>& y     \
        ) const                                                                 \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };

AtOp(min, (x.first() <= y.first() ? x : y))
AtOp(max, (x.first() >= y.first() ? x : y))
AtOp(minMagSqr, (magSqr(x.first())<=magSqr(y.first()) ? x : y))
AtOp(maxMagSqr, (magSqr(x.first())>=magSqr(y.first()) ? x : y))

#undef AtOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Three-way comparison operation of two parameters,
//  similar to the \c <=> operator in C++20.
//
//  \return a negative value for less, a positive value for greater,
//      and zero for equal value.
template<class T>
struct compareOp
{
    [[nodiscard]] __attribute__((always_inline))
    int  operator()(const T& a, const T& b) const
    {
        return (a < b) ? -1 : (b < a) ? 1 : 0;
    }
};

/**
 * An operator that applies multiple independent operators to the members of a tuple-esque-thing.
 *
 * For example, `ParallelOp<minOp<scalar>, sumOp<scalar>>` is an operator that produces the same
 * effect as applying `minOp` to `x[0]` and applying `sumOp` to `x[1]` (for a given input list
 * `x` of size >= 2).
 *
 * This is useful when you want to do collective MPI operations more efficiently. Instead of
 * something like this (which does 3 reductions in a row, repeating the whole communication process
 * three times):
 *
 *    reduce(data[0], minOp<scalar>());
 *    reduce(data[1], maxOp<scalar>());
 *    reduce(data[2], sumOp<scalar>());
 *
 * You can do:
 *    reduce(data, ParallelOp<minOp<scalar>, maxOp<scalar>, sumOp<scalar>>{});
 *
 *
 * Instead of something like:
 *
 *    scalar x, y = whatever...;
 *    reduce(x, sumOp<scalar>{});
 *    reduce(y, sumOp<scalar>{});
 *
 * You can do:
 *    reduce(std::tie(x, y), UniformParallelOp<sumOp<scalar>, 2>{});
 *
 * This version will do one round of communications (in which the messages are larger) to achieve the
 * same effect. Since most networking/etc. communication systems have the property that the constant
 * part of the cost dominates the linear part for small messages (for fairly large values of "small"),
 * the expectation is that this is a nice little win. Reduce operations are an extreme case of this, since
 * the thing being reduced is typically *a number*, so the cost of sending 3 of them is plausibly *equal*
 * to the cost of sending one.
 *
 * Note that this implementation will only work for stateless sub-operators, since it reconstructs them
 * at the point of use each time (to facilitate inlining of trivial operators). If your operator is a nontrivial
 * one that stores state in the operator struct, some adaptation of this will be required (just construct and
 * store the state in a field at construction-time). Since I can't find any such operators in the codebase -
 * despite the operator system being seemingly designed to handle them - I'm not bothering to resolve this issue.
 */
template<typename... Ops>
struct ParallelOp {
    // Special version for tuples.
    template<typename... Ts, size_t... Is>
    auto run(const std::tuple<Ts...>& x, const std::tuple<Ts...>& y, const std::index_sequence<Is...>) const {
        auto out = x;
        // Fun fact: `operator,` is a thing, and we can do foldexprs on it:
        ((std::get<Is>(out) = Ops{}(std::get<Is>(x), std::get<Is>(y))), ...);
        return out;
    };

    template<typename T, size_t... Is>
    auto run(const T& x, const T& y, const std::index_sequence<Is...>) const {
        T out = x;
        ((out[Is] = Ops{}(x[Is], y[Is])), ...);
        return out;
    };

    // If we had concepts we could actually constrain this a bit, but since we don't
    // we just gotta leave all the type variables free and hope the substitution errors
    // aren't too horrifying if someone uses this wrong by accident :D.
    template<typename T>
    auto operator()(const T& x, const T& y) const {
        return run(x, y, std::index_sequence_for<Ops...>{});
    }

    // Extract the I'th operator.
    template<size_t I>
    using getOp = std::tuple_element<I, std::tuple<Ops...>>;
};

// A fused operator is commutative if all of its member operators are commutative (since they
// run "in parallel").
template<typename... Ts>
constexpr bool IsCommutative<ParallelOp<Ts...>> = (IsCommutative<Ts> && ...);

template<size_t I, typename T>
using DiscardIndex = T;

// Non-executed function used to calculate the type for UniformParallelOp
template<typename Op, size_t... Is>
constexpr auto uniformParallelOpImpl(const std::index_sequence<Is...>&) {
    return ParallelOp<DiscardIndex<Is, Op>...>{};
}

/**
 * Same as ParallelOp<Op, Op, Op, ... (N times)>
 *
 * Handy for when you want to do the same thing to all the members of a list of constexpr length.
 */
template<typename Op, int N>
using UniformParallelOp = decltype(uniformParallelOpImpl<Op>(std::make_index_sequence<N>{}));

/// Determine if a ParallelOp<stuff...> is uniform (ie. all the stuff... are
/// the same).
template<typename... Ops>
constexpr bool isUniformParallelOp(const ParallelOp<Ops...>&) {
    if constexpr (sizeof...(Ops) < 2) {
        return true;
    } else {
        // Yuck.
        using First = std::tuple_element_t<0, std::tuple<Ops...>>;
        return (std::is_same_v<First, Ops> && ...);
    }
};
/// Martian operators are not uniform parallel ops.
template<typename T>
constexpr bool isUniformParallelOp(const T&) {
    return false;
}

template<typename>
struct isParallelOp: std::false_type {};

template<typename ...T>
struct isParallelOp<ParallelOp<T...>>: std::true_type {};

// Non-executed function.
template<typename... Ops>
constexpr auto getFirstOpImpl(const ParallelOp<Ops...>&) {
    return std::tuple_element_t<0, std::tuple<Ops...>>{};
}

/**
 * Get the first Op from a ParallelOp<...>.
 */
template<typename T>
using getFirstOp = decltype(getFirstOpImpl(T{}));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
