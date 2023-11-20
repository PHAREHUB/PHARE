#include "mpi.h"
#include <fstream>
#include <iostream>
#include <exception>

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

    std::ofstream file("PHARE_MPI_Get_library_version.txt", std::ios::out);
    if (!file)
        throw std::runtime_error("error opening version file");
    file << version << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    try
    {
        config_mpi();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return 0;
}
