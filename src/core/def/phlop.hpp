#ifndef PHARE_CORE_DEF_PHLOP_HPP
#define PHARE_CORE_DEF_PHLOP_HPP

#if __has_include("phlop/timing/mpi_scope_timer.hpp")

#if __has_include("mpi.h")
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
