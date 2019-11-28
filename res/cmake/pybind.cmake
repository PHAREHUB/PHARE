find_package(pybind11 CONFIG QUIET)

if (NOT pybind11_FOUND)

  message("downloading subproject pybind11")
  set(PYBIND11_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pybind11)

  if (NOT EXISTS ${PYBIND11_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/pybind/pybind11 ${PYBIND11_SRCDIR}
    )
  endif()

  add_subdirectory(${PYBIND11_SRCDIR})

endif()
