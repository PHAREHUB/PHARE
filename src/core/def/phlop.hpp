#ifndef PHARE_CORE_DEF_PHLOP_HPP
#define PHARE_CORE_DEF_PHLOP_HPP

#if __has_include("phlop/timing/mpi_scope_timer.hpp")

// The timer is picked per translation unit: phlop::mpi:: (MPI_Wtime) if mpi.h is on the
//  include path (amr, simulator, python), else serial phlop:: (core and core tests, no MPI).
//  They are separate instances, samrai.cpp only initialises/writes phlop::mpi::, so scopes
//  in TUs without MPI are not recorded there.
#if __has_include("mpi.h")
// phlop includes mpi.h itself, phare_mpi.hpp must come first: with -Dphare_configurator=ON
//  its generated core/def/_gen_mpi.hpp is what defines OMPI/MPICH_SKIP_MPICXX, without it
//  the MPI C++ bindings get compiled in and fail to link (no libmpi_cxx)
#include "phare_mpi.hpp"                    // IWYU pragma: keep
#include "phlop/timing/mpi_scope_timer.hpp" // IWYU pragma: keep
#else
#include "phlop/timing/scope_timer.hpp" // IWYU pragma: keep
#endif

#define PHARE_HAVE_PHLOP 1
#define PHARE_WITH_PHLOP(...) __VA_ARGS__

#define PHARE_SCOPE_TIMER PHLOP_SCOPE_TIMER
#define PHARE_SCOPE(str, ...) PHLOP_SCOPE(str, __VA_ARGS__)

// phlop scope options per log level, any of time, mem, perf_light, perf_all,
//  eg -DPHARE_LOG_SCOPE_FLAGS_1=time,mem,perf_light
//  mem and perf_* are only recorded when enabled at runtime:
//    PHLOP_SCOPE_MEM=1, PHLOP_PERF_LIGHT / PHLOP_PERF_ALL, see phlop/timing/scope.hpp
#ifndef PHARE_LOG_SCOPE_FLAGS_1
#define PHARE_LOG_SCOPE_FLAGS_1 PHLOP_SCOPE_DEFAULT
#endif
#ifndef PHARE_LOG_SCOPE_FLAGS_2
#define PHARE_LOG_SCOPE_FLAGS_2 PHLOP_SCOPE_DEFAULT
#endif
#ifndef PHARE_LOG_SCOPE_FLAGS_3 // hot loops, a memory sample costs ~10us
#define PHARE_LOG_SCOPE_FLAGS_3 time
#endif

#else
#define PHARE_HAVE_PHLOP 0
#define PHARE_WITH_PHLOP(...)

#endif // __has_include(...)

#endif /* PHARE_CORE_DEF_PHLOP_HPP */
