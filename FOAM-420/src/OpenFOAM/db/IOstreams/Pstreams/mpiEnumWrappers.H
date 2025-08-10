#ifndef FOAMCORE_MPIENUMWRAPPERS_H
#define FOAMCORE_MPIENUMWRAPPERS_H


namespace Foam {

// We need to metaprogram in terms of various MPI enumerants.
// Obnoxiously, many mpi implementations use non-constexpr-compatible
// types to represent datatypes and operators. So: we define our own versions here, use
// these to metaprogram, and then the 1-1 mapping will optimise out later :D.
// FUN.
enum class MpiOp {
    SUM,
    MIN,
    MAX,
    PROD,
    LOR,
    LAND,
    REPLACE
};

template<MpiOp op>
MPI_Op toMpiOp() {
    switch (op) {
        case MpiOp::SUM:
            return MPI_SUM;
        case MpiOp::MIN:
            return MPI_MIN;
        case MpiOp::MAX:
            return MPI_MAX;
        case MpiOp::PROD:
            return MPI_PROD;
        case MpiOp::LOR:
            return MPI_LOR;
        case MpiOp::LAND:
            return MPI_LAND;
        case MpiOp::REPLACE:
            return MPI_REPLACE;
    }
}

enum class MpiType {
    FLOAT,
    DOUBLE,
    CHAR,
    UNSIGNED_CHAR,
    INT,
    UNSIGNED,
    LONG,
    UNSIGNED_LONG,
    SHORT,
    UNSIGNED_SHORT,
    LONG_LONG,
    UNSIGNED_LONG_LONG,
    BYTE,
    UINT8_T,
    UINT16_T,
    UINT32_T,
    UINT64_T,
    INT8_T,
    INT16_T,
    INT32_T,
    INT64_T
};


template<MpiType t>
MPI_Datatype toMpiType() {
    switch (t) {
        case MpiType::FLOAT:
            return MPI_FLOAT;
        case MpiType::DOUBLE:
            return MPI_DOUBLE;
        case MpiType::CHAR:
            return MPI_CHAR;
        case MpiType::UNSIGNED_CHAR:
            return MPI_UNSIGNED_CHAR;
        case MpiType::INT:
            return MPI_INT;
        case MpiType::UNSIGNED:
            return MPI_UNSIGNED;
        case MpiType::LONG:
            return MPI_LONG;
        case MpiType::UNSIGNED_LONG:
            return MPI_UNSIGNED_LONG;
        case MpiType::SHORT:
            return MPI_SHORT;
        case MpiType::UNSIGNED_SHORT:
            return MPI_UNSIGNED_SHORT;
        case MpiType::LONG_LONG:
            return MPI_LONG_LONG;
        case MpiType::UNSIGNED_LONG_LONG:
            return MPI_UNSIGNED_LONG_LONG;
        case MpiType::BYTE:
            return MPI_BYTE;
        case MpiType::UINT8_T:
            return MPI_UINT8_T;
        case MpiType::UINT16_T:
            return MPI_UINT16_T;
        case MpiType::UINT32_T:
            return MPI_UINT32_T;
        case MpiType::UINT64_T:
            return MPI_UINT64_T;
        case MpiType::INT8_T:
            return MPI_INT8_T;
        case MpiType::INT16_T:
            return MPI_INT16_T;
        case MpiType::INT32_T:
            return MPI_INT32_T;
        case MpiType::INT64_T:
            return MPI_INT64_T;
    }
}


} // namespace Foam


#endif //FOAMCORE_MPIENUMWRAPPERS_H
