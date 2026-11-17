
# xsimd (header only): explicit SIMD backend of the hybrid particle push and deposit
#  -DwithXsimd=OFF to build without it, the lane loops are then vectorised by the compiler
#  the xsimd kernels are compiled only when the target has AVX-512 (e.g. -march=native on
#  a CPU that has it), see src/core/numerics/pusher/hyb_xsimd.hpp

set (PHARE_HAS_XSIMD "0")
if (withXsimd)

  if(NOT DEFINED PHARE_XSIMD_VERSION)
    SET(PHARE_XSIMD_VERSION "13.0.0")
  endif()

  set(XSIMD_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/xsimd)
  phare_github_get_or_update(xsimd ${XSIMD_SRCDIR} xtensor-stack/xsimd ${PHARE_XSIMD_VERSION})
  include_directories(SYSTEM ${XSIMD_SRCDIR}/include)

  set (PHARE_HAS_XSIMD "1")
endif(withXsimd)

add_compile_definitions(PHARE_WITH_XSIMD=${PHARE_HAS_XSIMD})
