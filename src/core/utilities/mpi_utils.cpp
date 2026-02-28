//

#include "core/errors.hpp"

#include "mpi_utils.hpp"

namespace PHARE::core::mpi
{
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


bool any(bool b)
{
    int global_sum, local_sum = static_cast<int>(b);
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return global_sum > 0;
}

bool any_errors()
{
    PHARE_DEBUG_DO({ return any(core::Errors::instance().any()); }) // ALLREDUCE!
    return false;
}


void log_error(std::string const key, std::string const val)
{
#ifdef NDEBUG
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    core::Errors::instance().log(key, val);
}

void barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}


void debug_barrier()
{
#if PHARE_DEBUG_BARRIER
    barrier();
#endif
}



std::string date_time(std::string format)
{
    std::time_t t = std::time(NULL);
    char buffer[80];
    struct tm ti;
    localtime_r(&t, &ti);
    std::strftime(buffer, 80, format.c_str(), &ti);
    std::string date_time{buffer};
    return collect(date_time)[0];
}

std::int64_t unix_timestamp_now()
{ // static cast as std::time_t typedef can vary across operating systems
    return all_get_from_rank_0([]() { return static_cast<std::int64_t>(std::time(NULL)); });
}

} // namespace PHARE::core::mpi
