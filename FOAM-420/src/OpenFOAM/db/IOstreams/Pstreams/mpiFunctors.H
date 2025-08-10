#ifndef MPI_FUNCTORS_H
#define MPI_FUNCTORS_H

// Trampoline functions (as functors) for calling various MPI functions.
//
// Makes it possible to pass equivalently-shaped MPI primitives as template parmaeters, which can reduce
// duplication a lot. Also serves as a layer of indirection between this library and MPI to babysit the insanity
// of the dynamic-code mechanism (which does not have reasonable include management).
//
// These would be better if defined as templates, but the "dynamic code" system's inability to handle mpi.h prevents this.
// Were that done, we could get some improved type-safety and TBAA properties out of this, and a teeny-tiny speedup:
// this version creates at least one redundant stack frame per call.


namespace Foam {

struct mpiAllGatherv {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, const int recvcounts[],
        const int displs[],
        int comm
    );
};

struct mpiAllGather {
    int operator()(
        const void* sendbuf, int  sendcount,
        void* recvbuf, int recvcount,
        int comm
    );
};

struct mpiGatherv {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, const int recvcounts[],
        const int displs[],
        int comm
    );
};

struct mpiGather {
    int operator()(
        const void* sendbuf, int sendcount,
        void* recvbuf, int recvcount,
        int comm
    );
};

} // Namespace Foam

#endif // MPI_FUNCTORS_H
