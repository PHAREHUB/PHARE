#include "mpi_utils.h"

namespace PHARE::core::mpi
{
bool initialized()
{
    int i;
    MPI_Initialized(&i);
    return i > 0;
}


int size()
{
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    return mpi_size;
}

int rank()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    return mpi_rank;
}

std::size_t max(std::size_t const local, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto perMPI = collect(local, mpi_size);
    return *std::max_element(std::begin(perMPI), std::end(perMPI));
}


void abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}

void finalize()
{
    MPI_Finalize();
}

bool any(bool b)
{
    int global_sum, local_sum = static_cast<int>(b);
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return global_sum > 0;
}

} // namespace PHARE::core::mpi
