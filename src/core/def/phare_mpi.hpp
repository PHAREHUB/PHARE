#ifndef PHARE_CORE_DEF_MPI_HPP
#define PHARE_CORE_DEF_MPI_HPP


// DO NOT INCLUDE MPI MANUALLY! USE THIS FILE!

#if __has_include("core/def/_gen_mpi.hpp")
#include "core/def/_gen_mpi.hpp"
#else
// Not always an issue, but not recommended
#endif

#ifndef OMPI_SKIP_MPICXX
// avoids default including mpicxx
#define OMPI_SKIP_MPICXX 1
#endif /* OMPI_SKIP_MPICXX */

#ifndef MPICH_SKIP_MPICXX
// avoids default including mpicxx
#define MPICH_SKIP_MPICXX 1
#endif /* MPICH_SKIP_MPICXX */

#include "core/def/pragma_disable.hpp"

// clang-format off
DISABLE_WARNING(cast-function-type, bad-function-cast, 42)
#include "mpi.h"
ENABLE_WARNING(cast-function-type, bad-function-cast, 42)
// clang-format on


#endif /*PHARE_CORE_DEF_MPI_HPP*/
