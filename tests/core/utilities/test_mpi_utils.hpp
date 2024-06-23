#ifndef PHARE_TEST_CORE_UTILITIES_MPI_HPP
#define PHARE_TEST_CORE_UTILITIES_MPI_HPP

#include "core/utilities/mpi_utils.hpp"
#include <iostream>

namespace PHARE::core::mpi
{


struct Lifecycle
{
    static inline bool error_occurred = false;

    Lifecycle(int argc, char** argv)
    {
        int err = MPI_Init(&argc, &argv);
        if (err != MPI_SUCCESS)
        {
            std::cerr << "MPI Initialization failed with error code: " << err << std::endl;
            error_occurred = true;
        }
    }
    ~Lifecycle()
    {
        if (!error_occurred)
        {
            int err = MPI_Finalize();
            if (err != MPI_SUCCESS)
            {
                std::cerr << "MPI Finalization failed with error code: " << err << std::endl;
            }
        }
    }
};

} // namespace PHARE::core::mpi


#endif /* PHARE_TEST_CORE_UTILITIES_MPI_H */
