
#ifndef PHARE_CORE_DEF_KOKKOS_TOOLS_HPP
#define PHARE_CORE_DEF_KOKKOS_TOOLS_HPP

#if __has_include(<Kokkos_Core.hpp>)
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#define PHARE_HAVE_KOKKOS 1
#define PHARE_WITH_KOKKOS(...) __VA_ARGS__

#else

#define PHARE_HAVE_KOKKOS 0
#define PHARE_WITH_KOKKOS(...)

#endif //


#if PHARE_HAVE_KOKKOS && __has_include("KokkosTools/profiling/all/kp_all.hpp")

#include "KokkosTools/profiling/all/kp_all.hpp"

#define PHARE_HAVE_KOKKOS_TOOLS 1
#define PHARE_WITH_KOKKOS_TOOLS(...) __VA_ARGS__
#else
#define PHARE_HAVE_KOKKOS_TOOLS 0
#define PHARE_WITH_KOKKOS_TOOLS(...)

#endif // __has_include(...)

#endif /* PHARE_CORE_DEF_KOKKOS_TOOLS_HPP */
