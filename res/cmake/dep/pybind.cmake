
SET(PYBIND_MIN_VERSION "2.5.0")

function(get_pybind)

  # Pybind errors with clang, it is default in GCC
  # https://github.com/pybind/pybind11/issues/1604
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (PHARE_FLAGS ${PHARE_FLAGS} -fsized-deallocation)
  endif()

  message("downloading subproject pybind11")
  set(PYBIND11_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pybind11)
  set(PYBIND11_VERSION master)

  phare_github_get_or_update(pybind11 ${PYBIND11_SRCDIR} pybind/pybind11 ${PYBIND11_VERSION})

  add_subdirectory(${PYBIND11_SRCDIR})

endfunction(get_pybind)

if (forceGetPybind)
  get_pybind()
else()

  find_package(pybind11 ${PYBIND_MIN_VERSION} CONFIG QUIET)

  if (NOT pybind11_FOUND)
    get_pybind()
  endif()

endif(forceGetPybind)
