#include "mpi_utils.h"

namespace PHARE::core::mpi {



std::vector<std::string> collectStrings(std::string str, int mpi_size = 0,
                                        std::string null_str = "null")
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    str             = str.empty() ? null_str : str;
    auto maxMPISize = max(str.size(), mpi_size);
    auto perMPI     = collect(str.size(), mpi_size);
    std::vector<char> chars(maxMPISize * mpi_size);
    _collect(str.c_str(), chars, str.size(), maxMPISize);
    std::vector<std::string> values;
    for (int i = 0; i < mpi_size; i++)
    {
        std::string data{&chars[maxMPISize * i], perMPI[i]};
        data = data == null_str ? "" : data;
        values.emplace_back(data);
    }
    return values;
}



size_t max(size_t local, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    auto perMPI = collect(local, mpi_size);
    return *std::max_element(std::begin(perMPI), std::end(perMPI));
}
}
