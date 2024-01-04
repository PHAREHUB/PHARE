#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1

#include "mpi.h"
#include <fstream>
#include <iostream>
#include <exception>
#include <string_view>

#if __has_include("hdf5.h")
#include "hdf5.h"
constexpr static std::string_view hdf5_version = H5_VERSION;
#else
constexpr static std::string_view hdf5_version     = "HDF5 was not found!";
constexpr static std::string_view H5_HAVE_PARALLEL = false;
#endif
constexpr std::string_view hdf5_is_parallel()
{
    if constexpr (H5_HAVE_PARALLEL)
        return "yes";
    return "no";
}


void write_string_to_file(std::string const& buff, std::string const& filename)
{
    std::ofstream file(filename, std::ios::out);
    if (!file)
        throw std::runtime_error(std::string{"error opening file: "} + filename);
    file << buff << std::endl;
}

void write_hdf5_version()
{
    write_string_to_file(std::string{hdf5_version}, "PHARE_HDF5_version.txt");
}

void write_hdf5_is_parallel()
{
    write_string_to_file(std::string{hdf5_is_parallel()}, "PHARE_HDF5_is_parallel.txt");
}

void write_hdf5_info()
{
    write_hdf5_version();
    write_hdf5_is_parallel();
}

void config_mpi()
{
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        throw std::runtime_error("MPI not initialized");
    }

    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int resultlen = 0;
    int ret       = MPI_Get_library_version(version, &resultlen);
    if (ret != MPI_SUCCESS)
        throw std::runtime_error("error calling MPI_Get_library_version");

    write_string_to_file(version, "PHARE_MPI_Get_library_version.txt");
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    try
    {
        config_mpi();
        write_hdf5_info();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return 0;
}
