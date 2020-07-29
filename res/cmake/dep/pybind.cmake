

function(get_pybind)

  # Pybind errors with clang, it is default in GCC
  # https://github.com/pybind/pybind11/issues/1604
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (PHARE_FLAGS ${PHARE_FLAGS} -fsized-deallocation)
  endif()

  message("downloading subproject pybind11")
  set(PYBIND11_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pybind11)

  if (NOT EXISTS ${PYBIND11_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/pybind/pybind11 ${PYBIND11_SRCDIR} --depth 1 -b master
    )
  else()
    if(devMode)
      message("downloading latest pybind11 updates")
      execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${PYBIND11_SRCDIR})
    endif(devMode)
  endif()

  add_subdirectory(${PYBIND11_SRCDIR})

endfunction(get_pybind)

if (forceGetPybind)
  get_pybind()
else()

  find_package(pybind11 CONFIG QUIET)

  if (NOT pybind11_FOUND)
    get_pybind()
  endif()

endif(forceGetPybind)
