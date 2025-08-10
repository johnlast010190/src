#include <mpi.h>

#include "mpiFunctors.H"
#include "PstreamGlobals.H"

// This project has a "fun" convention of randomly including `.C` files into headers to improve
// inter-procedural optimisation (and presumably to warm people's offices). Ordinarily such include
// files would have a `.hpp` or `.inc` extension, but oh well.
//
// Doing that with this file will result in Unspecified Bad Things.

bool buffersDoNotAlias(
    const void* x, int xLen,
    const void* y, int yLen
) {
    uintptr_t X = (uintptr_t) x;
    uintptr_t Y = (uintptr_t) y;

    return (X >= (Y + yLen)) || ((X + xLen) <= Y);
}

int Foam::mpiAllGatherv::operator()(
    const void* sendbuf, int sendcount,
    void* recvbuf, const int recvcounts[],
    const int displs[], int comm
) {
    return MPI_Allgatherv(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcounts, displs, MPI_BYTE, PstreamGlobals::MPICommunicators_[comm]);
}

int Foam::mpiAllGather::operator()(
    const void* sendbuf, int sendcount,
    void* recvbuf, int recvcount,
    int comm
) {
    if (sendbuf == nullptr) {
        sendbuf = MPI_IN_PLACE;
    }

    FOAM_ASSERT(sendbuf == MPI_IN_PLACE || buffersDoNotAlias(sendbuf, sendcount, recvbuf, recvcount)) {
        FatalErrorInFunction << "MPI allgather buffers may not alias" << Foam::abort(FatalError);
    }

    return MPI_Allgather(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcount, MPI_BYTE, PstreamGlobals::MPICommunicators_[comm]);
}

int Foam::mpiGatherv::operator()(
    const void* sendbuf, int sendcount,
    void* recvbuf, const int recvcounts[],
    const int displs[], int comm
) {
    return MPI_Gatherv(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcounts, displs, MPI_BYTE, 0, PstreamGlobals::MPICommunicators_[comm]);
}

int Foam::mpiGather::operator()(
    const void* sendbuf, int sendcount,
    void* recvbuf, int recvcount, int comm
) {
    if (sendbuf == nullptr) {
        sendbuf = MPI_IN_PLACE;
    }

    FOAM_ASSERT(sendbuf == MPI_IN_PLACE || buffersDoNotAlias(sendbuf, sendcount, recvbuf, recvcount)) {
        FatalErrorInFunction << "MPI gather buffers may not alias" << Foam::abort(FatalError);
    }

    return MPI_Gather(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcount, MPI_BYTE, 0, PstreamGlobals::MPICommunicators_[comm]);
}
