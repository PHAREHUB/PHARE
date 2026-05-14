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
#else
#define PHARE_HAVE_PHLOP 0
#define PHARE_WITH_PHLOP(...)

#endif // __has_include(...)

#endif /* PHARE_CORE_DEF_PHLOP_HPP */
