#ifndef FOAMCORE_MPIRequest_H
#define FOAMCORE_MPIRequest_H

#include <mpi.h>


/**
 * Represents a result value that will eventually exist for an MPI request.
 * @tparam T
 */
class MPIRequest {
    MPI_Request req{MPI_REQUESI_NULL};

public:

    /// Wait for the underlying request to finish, and return the result value.
    /// This may not be called from multiple threads, although that could be changed.
    /// This may not be called twice.
    void await() {
        if (req == MPI_REQUESI_NULL) {
            return;
        } else {
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }
    }

    /**
     * Checks if the request has finished yet or not. Allows you to test if the value
     * is opportunistically ready to read without having to block.
     *
     * Use of this function is discouraged, since it is inherently a race condition. It is
     * fairly easy to accidentally introduce nondeterminism into your program if you use
     * this function. It is also possible to make things faster, if done correctly :D.
     */
    bool isReady() {
        int out;
        MPI_Request_get_status(req, &out, MPI_STATUS_IGNORE);
        return out;
    }

    operator MPI_Request*() {
        return &req;
    }
};



#endif //FOAMCORE_MPIRequest_H
