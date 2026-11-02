
if (withKokkos)

  set(KOKKOS_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/Kokkos)
  set(KOKKOS_TOOLS_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/KokkosTools)

  phare_github_get_or_update(Kokkos ${KOKKOS_SRCDIR} kokkos/kokkos develop)
  phare_github_get_or_update(KokkosTools ${KOKKOS_TOOLS_SRCDIR} kokkos/kokkos-tools develop)

  set(BUILD_SHARED_LIBS ON CACHE BOOL "Make shared libs" FORCE)

  add_subdirectory(${KOKKOS_SRCDIR})
  include_directories(${KOKKOS_SRCDIR}/core/src)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/subprojects/Kokkos)

  add_subdirectory(${KOKKOS_TOOLS_SRCDIR})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/subprojects/KokkosTools/common)

  set(PHARE_BASE_LIBS "${PHARE_BASE_LIBS}" Kokkos::kokkoscore)

endif(withKokkos)
