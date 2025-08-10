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
    (c) 2016 OpenCFD Ltd.

InNamespace
    Foam

Description
    Inter-processor communication reduction functions.

\*---------------------------------------------------------------------------*/

#ifndef PstreamReduceOps_H
#define PstreamReduceOps_H

#include <mpi.h>
#include "primitives/ops/ops.H"
#include "primitives/Vector2D/vector2D/vector2D.H"
#include "PstreamGlobals.H"
#include "UPstream.H"
#include "Pstream.H"
#include "MPIRequest.H"
#include "mpiEnumWrappers.H"
#include "primitives/Vector/Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename>
struct is_tuple: std::false_type {};

template<typename ...T>
struct is_tuple<std::tuple<T...>>: std::true_type {};


/**
 * Make a normal binary operator template usable by MPI's horrible C API.
 *
 * While an instance of this class template is in scope, the MPI_Op value exposed by
 * getOp() can be passed to MPI functions such as `reduce()`. A reasonable way to
 * use this class is to just construct an instance on the stack right before you call
 * some MPI operation, such that this gets neatly cleaned up right after. You could
 * save a teensy bit of overhead (approximately 2 heap allocations) if you hoist
 * the scope of the operator declaration and reuse it.
 *
 * Obviously, make sure you don't deallocate the object until you're done :D
 */
template<typename Op, typename T>
class CustomMPIReduceOp {
    MPI_Op op;

public:
    MPI_Op getOp() {
        return op;
    }

    CustomMPIReduceOp() {
        MPI_Op_create(&CustomMPIReduceOp<Op, T>::doTheThing, IsCommutative<Op>, &op);
    }
    ~CustomMPIReduceOp() {
        MPI_Op_free(&op);
    }

    /**
     * Of course, MPI requires that we represent our function in the craziest possible way.
     *
     * Ignoring the datatype provided by MPI here is not an accident: we propagate the type information ourselves
     * using C++ types (in the template parameter), and as such we guarantee that an instantiation of this template
     * with a type T is only used by an MPI call using a compatible type, so we can safely assume that the `void*`
     * pointers secretly point to a T.
     *
     * `len` being a pointer appears to just be a deranged quirk of how MPI defines its operator contract. Since
     * MPI's structure prevents this form ever inlining, this will always be accessed as a call through a function
     * pointer (which is quite inefficient), so this is mitigated by this hilarious structure where the operator
     * has to pairwise reduce two lists, instead of just being a trivial reducer that gets inlined everywhere.
     *
     * Sigh.
     */
    static void doTheThing(void* invec, void* inoutvec, int* len, MPI_Datatype*) {
        const T* in = (T*) invec;
        T* out = (T*) inoutvec;
        for (int i = 0; i < *len; i++) {
            out[i] = Op{}(in[i], out[i]);
        }
    }
};

template<bool ToMaster, bool Async, typename T>
MPIRequest performReduce
(
    T& Value,
    MPI_Op MPIOp,
    MPI_Datatype MPIType,
    int MPICount,
    const label comm
)
{
    MPIRequest out;
    if constexpr (ToMaster) {
        bool isMaster = UPstream::master(comm);

        if constexpr (Async) {
            MPI_IReduce(
                isMaster ? MPI_IN_PLACE : &Value,
                &Value,
                MPICount,
                MPIType,
                MPIOp,
                UPstream::masterNo(),
                PstreamGlobals::MPICommunicators_[comm],
                out
            );
        } else {
            MPI_Reduce(
                isMaster ? MPI_IN_PLACE : &Value,
                &Value,
                MPICount,
                MPIType,
                MPIOp,
                UPstream::masterNo(),
                PstreamGlobals::MPICommunicators_[comm]
            );
        }
    } else {
        if constexpr (Async) {
            MPI_Iallreduce(
                MPI_IN_PLACE,
                &Value,
                MPICount,
                MPIType,
                MPIOp,
                PstreamGlobals::MPICommunicators_[comm],
                out
            );
        } else {
            MPI_Allreduce(
                MPI_IN_PLACE,
                &Value,
                MPICount,
                MPIType,
                MPIOp,
                PstreamGlobals::MPICommunicators_[comm]
            );
        }
    }
    return out;
}

/// Convert a C++ builtin type to an MPI type.
/// This entire thing should be constexpr-if, but the OpenMPI people decided to put TBAA-unsafe casts
/// in their definitions of the type enums, meaning they cannot be used in a constexpr context (among
/// other things). Cry.
/// Still: the constant propagator should have fun with this function.
template<typename T>
__attribute__((always_inline))
constexpr MpiType getMpiTypeForBuiltinType() {
    if constexpr (std::is_same_v<float, T>) {
        return MpiType::FLOAT;
    } else if constexpr (std::is_same_v<double, T>) {
        return MpiType::DOUBLE;
    } else if constexpr (std::is_same_v<char, T>) {
        return MpiType::CHAR;
    }  else if constexpr (std::is_same_v<unsigned char, T>) {
        return MpiType::UNSIGNED_CHAR;
    } else if constexpr (std::is_same_v<int, T>) {
        // MPI specifies that its types refer to these variants of the C++ types, not the
        // `uint32_t` fixed-width ones. In the unlikely event the rest of FOAM is written to
        // work on standards-conforming C++ implementations with differing int sizes, this
        // pedantry may prove useful.
        return MpiType::INT;
    } else if constexpr (std::is_same_v<unsigned int, T>) {
        return MpiType::UNSIGNED;
    } else if constexpr (std::is_same_v<long, T>) {
        return MpiType::LONG;
    } else if constexpr (std::is_same_v<unsigned long, T>) {
        return MpiType::UNSIGNED_LONG;
    } else if constexpr (std::is_same_v<short, T>) {
        return MpiType::SHORT;
    } else if constexpr (std::is_same_v<unsigned short, T>) {
        return MpiType::UNSIGNED_SHORT;
    } else if constexpr (std::is_same_v<long long, T>) {
        return MpiType::LONG_LONG;
    } else if constexpr (std::is_same_v<unsigned long long, T>) {
        return MpiType::UNSIGNED_LONG_LONG;
    } else if constexpr (std::is_same_v<uint8_t, T>) {
        return MpiType::UINT8_T;
    }  else if constexpr (std::is_same_v<uint16_t, T>) {
        return MpiType::UINT16_T;
    } else if constexpr (std::is_same_v<uint32_t, T>) {
        return MpiType::UINT32_T;
    } else if constexpr (std::is_same_v<uint64_t, T>) {
        return MpiType::UINT64_T;
    } else if constexpr (std::is_same_v<int8_t, T>) {
        return MpiType::INT8_T;
    } else if constexpr (std::is_same_v<int16_t, T>) {
        return MpiType::INT16_T;
    } else if constexpr (std::is_same_v<int32_t, T>) {
        return MpiType::INT32_T;
    } else if constexpr (std::is_same_v<int64_t, T>) {
        return MpiType::INT64_T;
    } else {
        // This'll let us move T around, but not much else.
        return MpiType::BYTE;
    }
}

template<typename TupleT, size_t... Is>
constexpr bool tupleIsUniformImpl(const std::index_sequence<Is...>&) {
    using FirstT = std::tuple_element_t<0, TupleT>;
    return (std::is_same_v<std::tuple_element_t<Is, TupleT>, FirstT> && ...);
}

/// Determine if all the members of a tuple are the same type.
template<typename TupleT>
constexpr bool tupleIsUniform() {
    constexpr int TupleLen = std::tuple_size_v<TupleT>;
    if constexpr (TupleLen == 0) {
        return true; // Errrr
    } else {
        return tupleIsUniformImpl<TupleT>(std::make_index_sequence<TupleLen>());
    }
}

/**
 * T should be described to MPI as a <N>-vector of elements of what type.
 *
 * This allows us to make use of the MPI builtin reduce operators in the maximum number of cases, provided
 * we describe our types appropriately. Add special cases here for every type that can be correctly reduced
 * by regarding it as an N-tuple of some builtin MPI type (and you want N independent reductions - one per
 * element).
 *
 * Builtin MPI operators are in
 * herently a bit more efficient because they don't need to use the silly
 * function-pointer-call nonsense (and in some cases have more elaborate things).
 */
template<typename T>
__attribute__((always_inline))
constexpr std::tuple<MpiType, int> mpiTypeFor() {
    static_assert(contiguous<T>(), "The given datatype is not contiguous, and MPI doesn't generally know how to do custom serialisation");
    static_assert(sizeof(T) <= MaxMpiCustomTypeLengthInBytes, "That type is too long for an MPI reduction. Increase MaxMpiCustomTypeLengthInBytes.");
    MpiType builtinT = getMpiTypeForBuiltinType<T>();

    // Spectacularly, MPI's enums rely on void* casts, so these cannot be constexpr-if. We trust, however, that
    // the constant propagator will save us...
    if (builtinT != MpiType::BYTE) {
        // T is built-in type, no further thought required.
        return {builtinT, 1};
    } else if (std::is_same_v<Vector2D<float>, T>) {
        return {MpiType::FLOAT, 2};
    } else if (std::is_same_v<Vector2D<double>, T>) {
        return {MpiType::DOUBLE, 2};
    } else {
        if constexpr (is_tuple<T>::value) {
            // Two nested constexpr-ifs because `tupleIsUniform` doesn't compile for things
            // that aren't tuples, and the inner one prevents the recursive call existing
            // after the second level (so the inliner doesn't get stuck).
            if constexpr (tupleIsUniform<T>()) {
                // A tuple where everything is the same type, and that one type is an MPI type,
                // can be handled by MPI quite nicely.
                auto eltT = mpiTypeFor<std::tuple_element_t<0, T>>();
                MpiType theT = std::get<0>(eltT);
                if (theT != MpiType::BYTE) {
                    return {theT, std::tuple_size_v<T>};
                }
            }
        }

        // No builtin representation exists.
        return {MpiType::BYTE, 0};
    }
}

/**
 * What is the builtin MPI operator for applying `Op` to elements of type `T`?
 *
 * Returns MPI_REPLACE if no such builtin exists.
 */
template<typename BinaryOp, typename T>
__attribute__((always_inline))
constexpr MpiOp getMpiOpForFunctor() {
    // The following builtin MPI ops are not used: BAND, BOR, LXOR, BXOR, MAXLOC, MINLOC, REPLACE, NOP
    if constexpr (isUniformParallelOp(BinaryOp{})) {
        // If it's a uniform parallel op, recurse into it and see if we can use MPI's builtin support
        // for multiple parallel reductions.
        return getMpiOpForFunctor<getFirstOp<BinaryOp>, T>();
    } else if constexpr (is_sumOp<BinaryOp>::value) {
        return MpiOp::SUM;
    } else if constexpr (is_minOp<BinaryOp>::value) {
        return MpiOp::MIN;
    } else if constexpr (is_maxOp<BinaryOp>::value) {
        return MpiOp::MAX;
    } else if constexpr (is_multiplyOp<BinaryOp>::value) {
        return MpiOp::PROD;
    } else if constexpr (is_orOp<BinaryOp>::value) {
        return MpiOp::LOR;
    } else if constexpr (is_andOp<BinaryOp>::value) {
        return MpiOp::LAND;
    } else {
        return MpiOp::REPLACE; // Used to represent it's a custom operation of some sort.
    }
}

template<typename BinaryOp, typename T>
void validateReduction(const T& v);

template<typename... Ts, typename... Ops, size_t... Is>
void validateTupleReduction(const std::tuple<Ts...>& Values, ParallelOp<Ops...>&&, std::index_sequence<Is...>&&) {
    // A parallel tuple reduction is valid if the equivalent set of independent non-batched reductions are
    // valid. This calls `validateReduction()` for each of those hypothetical reductions, to check them all
    // in turn.
    return (validateReduction<typename ParallelOp<Ops...>::template getOp<Is>>(std::get<Is>(Values)), ...);
}

// Abort compilation if reducing T with operator BinaryOp doesn't make sense as a thing to do.
// MPI isn't type-safe, so it'll just happily mangle the bytes for us, so we need to implement this
// compile-time check ourselves.
template<typename BinaryOp, typename T>
void validateReduction(const T& v) {
    if constexpr (is_tuple<T>::value) {
        if constexpr (isParallelOp<BinaryOp>::value) {
            return validateTupleReduction(v, BinaryOp{}, std::make_index_sequence<std::tuple_size_v<T>>{});
        }

        // For basic operators, you should use `op<T>` to reduce `T`.
        // If you've written your own operator: hopefully you wrote it correctly :D.
    } else if constexpr (is_sumOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, sumOp<T>>, "Mismatching operator type for MPI operation.");
    } else if constexpr (is_minOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, minOp<T>>, "Mismatching operator type for MPI operation.");
    } else if constexpr (is_maxOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, maxOp<T>>, "Mismatching operator type for MPI operation.");
    } else if constexpr (is_multiplyOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, multiplyOp<T>>, "Mismatching operator type for MPI operation.");
    } else if constexpr (is_orOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, orOp<T>>, "Mismatching operator type for MPI operation.");
    } else if constexpr (is_andOp<BinaryOp>::value) {
        static_assert(std::is_same_v<BinaryOp, andOp<T>>, "Mismatching operator type for MPI operation.");
    }
}

template<bool ToMaster, bool Async, typename T, typename BinaryOp>
MPIRequest reduceImpl
(
    T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
)
{
    // Doesn't actually do anything at runtime :D.
    validateReduction<BinaryOp>(Value);

    // Reduction is a no-op when not operating in parallel.
    if (!UPstream::parRun() || UPstream::nProcs(comm) == 1) {
        // Mildly horrifying
        return MPIRequest{};
    }

    // Firstly: are we sticking crayons up our nose or not?
    // It turns out there's some mildly deranged code that performs reductions over complex
    // structs of dynamic size, containing list<list<list<scalar>>> and other such delights.
    // There's no reasonable way to do that, so we must fall back to the original implementation
    // which does all the communication manually with point-to-point calls until we get around to
    // rewriting that code.
    // Fortunately, this check happens at compile-time (yay templates), and the slow path is taken
    // only very rarely. It just means we can't - yet - entirely delete the insane scatter/gather
    // functions.
    // A useful exercise for finding slow code is to add a `static_assert(contiguous<T>())` here,
    // and to inspect all of the things that then refuse to compile: those will be the slowest
    // reductions, and if they can be rephrased in terms of saner reductions then a large win may
    // be possible there.
    if constexpr (!contiguous<T>()) {
        static_assert(!is_tuple<T>::value);
        const List<Pstream::commsStruct>& comms =
                UPstream::nProcs(comm) < UPstream::nProcsSimpleSum ?
                    UPstream::linearCommunication(comm) :
                    UPstream::treeCommunication(comm);

        Pstream::gather(comms, Value, bop, Pstream::msgType(), comm);

        if constexpr (!ToMaster) {
            // We might actually be able to replace this step with a sane BCast, however...
            Pstream::scatter(comms, Value, Pstream::msgType(), comm);
        }

        // Mildly horrifying
        return MPIRequest{};
    } else {
        constexpr auto t = mpiTypeFor<T>();
        constexpr MpiType TheType = std::get<0>(t);

        // Attempt to match the reduction to one of MPI's built-in reduction operators, since
        // these are slightly more efficient (no need for awful function pointer calls at
        // every step).
        constexpr MpiOp builtinOpEquivalent = getMpiOpForFunctor<BinaryOp, T>();

        // The MPI_BYTE type coming out of `mpiTypeFor` represents a custom type. So if we're
        // doing something like `minOp<Vector>`, we will go to the other branch (since we need
        // the custom operator logic in that event).
        // NO_OP is used as a sentinel for "Some random custom operation".
        // Note that if you want to reduce 8-bit quantities, use MPI_UINT8_T, MPI_CHAR, etc.
        if constexpr (builtinOpEquivalent != MpiOp::REPLACE && TheType != MpiType::BYTE) {
            int width = std::get<1>(t);
            MPI_Datatype type = toMpiType<TheType>();
            MPI_Op builtinMpiOp = toMpiOp<builtinOpEquivalent>();
            return performReduce<ToMaster, Async>(Value, builtinMpiOp, type, width, comm);
        } else {
            // RAII a custom op for the duration of the operation.
            CustomMPIReduceOp<BinaryOp, T> op;

            // Treat it as a single record of sizeof(T)-many bytes. This'll get MPI to move it around how
            // we want, but it's up to our custom operator to know how to actually do the calculation.
            return performReduce<ToMaster, Async>(Value, op.getOp(), UPstream::dataTypes[sizeof(T)], 1, comm);
        }
    }
}

template<typename T, typename BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
) {
    reduceImpl<false, false>(Value, bop, comm);
}

template<typename T, typename BinaryOp>
auto reduceAsync
(
    T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
) {
    return reduceImpl<false, true>(Value, bop, comm);
}

template<typename T, typename BinaryOp>
void reduceToMaster
(
    T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
) {
    reduceImpl<true, false>(Value, bop, comm);
}

template<typename T, typename BinaryOp>
auto reduceToMasterAsync
(
    T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
) {
    return reduceImpl<true, true>(Value, bop, comm);
}

// Variant of the above that does orthogonal reductions across all elements of a list.
// This is a separate function template to allow for doing reductions that treat lists as
// elements (with some fancy operator), which would otherwise be unrepresentable if this
// were a specialisation of the above.
template<typename T, typename BinaryOp>
void listReduce
(
    UList<T>& list,
    const BinaryOp&,
    const label comm = UPstream::worldComm
)
{
    // Reduction is a no-op when not operating in parallel.
    if (!UPstream::parRun() || UPstream::nProcs(comm) == 1) {
        return;
    }
    static_assert(contiguous<T>());

    constexpr auto t = mpiTypeFor<T>();
    constexpr MpiType TheType = std::get<0>(t);

    // Attempt to match the reduction to one of MPI's built-in reduction operators, since
    // these are slightly more efficient (no need for awful function pointer calls at
    // every step).
    // Once again, we pray to the constant propagator because the OpenMPI devs are lunatics.
    constexpr MpiOp builtinOpEquivalent = getMpiOpForFunctor<BinaryOp, T>();

    // The MPI_BYTE type coming out of `mpiTypeFor` represents a custom type. So if we're
    // doing something like `minOp<Vector>`, we will go to the other branch (since we need
    // the custom operator logic in that event).
    // NO_OP is used as a sentinel for "Some random custom operation".
    // Note that if you want to reduce 8-bit quantities, use MPI_UINT8_T, MPI_CHAR, etc.
    if constexpr (builtinOpEquivalent != MpiOp::REPLACE && TheType != MpiType::BYTE) {
        MPI_Op builtinMpiOp = toMpiOp<builtinOpEquivalent>();
        MPI_Datatype type = toMpiType<TheType>();
        int width = std::get<1>(t);

        performReduce<false, false>(list[0], builtinMpiOp, type, width * list.size(), comm);
    } else {
        // RAII a custom op for the duration of the operation.
        CustomMPIReduceOp<BinaryOp, T> op;

        // Treat it as a single record of sizeof(T)-many bytes. This'll get MPI to move it around how
        // we want, but it's up to our custom operator to know how to actually do the calculation.
        performReduce<false, false>(list[0], op.getOp(), UPstream::dataTypes[sizeof(T)], list.size(), comm);
    }
}

template<bool ToMaster, bool Async, typename BinaryOp, size_t... Is, typename... Args>
void reduceTupleImpl (
    std::tuple<Args&...> Values,
    const BinaryOp& bop,
    const label comm,
    const std::index_sequence<Is...>&
) {
    // Most of the underlying operators probably do not work correctly with references,
    // and most of the things we want to reduce are not heap-ful anyway, so copying them
    // here is probably fine.
    std::tuple<Args...> copies{std::get<Is>(Values)...};
    reduceImpl<ToMaster, Async>(copies, bop, comm);
    ((std::get<Is>(Values) = std::get<Is>(copies)), ...);
}

/**
 * Special version of `reduce` for reducing a tuple of references with a composite operator.
 *
 * This is a mechanism for doing multiple independent reductions simultaneously, without
 * having to actually build a list or tuple upfront. ie:
 *
 * ```
 * scalar x = ...
 * scalar y = ...
 *
 * // Bad: slow
 * reduce(x, minOp<scalar>);
 * reduce(y, maxOp<scalar>);
 *
 * // Equivalent, but faster:
 * reduce(std::tie(x, y), ParallelOp<minOp<scalar>, maxOp<scalar>>);
 * ```
 */
template<typename BinaryOp, typename... Args>
void reduce
(
    std::tuple<Args&...> Values,
    const BinaryOp&,
    const label comm = UPstream::worldComm
) {
    reduceTupleImpl<false, false>(Values, BinaryOp{}, comm, std::index_sequence_for<Args...>{});
}
template<typename BinaryOp, typename... Args>
void reduceAsync
(
    std::tuple<Args&...> Values,
    const BinaryOp&,
    const label comm = UPstream::worldComm
) {
    return reduceTupleImpl<false, true>(Values, BinaryOp{}, comm, std::index_sequence_for<Args...>{});
}

template<typename BinaryOp, typename... Args>
void reduceToMaster
(
    std::tuple<Args&...> Values,
    const BinaryOp&,
    const label comm = UPstream::worldComm
) {
    reduceTupleImpl<true, false>(Values, BinaryOp{}, comm, std::index_sequence_for<Args...>{});
}

template<typename BinaryOp, typename... Args>
void reduceToMasterAsync
(
    std::tuple<Args&...> Values,
    const BinaryOp&,
    const label comm = UPstream::worldComm
) {
    return reduceTupleImpl<true, true>(Values, BinaryOp{}, comm, std::index_sequence_for<Args...>{});
}

// Reduce using either linear or tree communication schedule
template<class T, class BinaryOp>
T returnReduce
(
    const T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
)
{
    T WorkValue(Value);
    reduce(WorkValue, bop, comm);
    return WorkValue;
}

template<class T, class BinaryOp>
T returnReduceToMaster
(
    const T& Value,
    const BinaryOp& bop,
    const label comm = UPstream::worldComm
)
{
    T WorkValue(Value);
    reduceToMaster(WorkValue, bop, comm);
    return WorkValue;
}


// Reduce with sum of both value and count (for averaging)
template<class T>
void sumReduce
(
    T& Value,
    label& Count,
    const label comm = UPstream::worldComm
)
{
    reduce(std::tie(Value, Count), ParallelOp<sumOp<T>, sumOp<label>>{}, comm);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
