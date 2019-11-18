

if (NOT DEFINED PHARE_MPI_PROCS)
  set(PHARE_MPI_PROCS 2) # default MPI processes
endif()

# Pybind errors with clang, it is default in GCC
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS -fsized-deallocation)
endif()

if(coverage AND NOT MSVC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -DHAVE_EXECINFO_H -g3")
endif()
if(testMPI)
  function(add_phare_test binary directory)
    add_test(NAME ${binary} COMMAND mpirun -n ${PHARE_MPI_PROCS} ./${binary} WORKING_DIRECTORY ${directory})
    set_tests_properties(${binary} PROPERTIES ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
  endfunction(add_phare_test)
else()
  function(add_phare_test binary directory)
    add_test(NAME ${binary} COMMAND ${binary} WORKING_DIRECTORY ${directory})
  endfunction(add_phare_test)
endif(testMPI)

if(MSVC)
  set (PHARE_WERROR_FLAGS /W4 /WX)
else()
  set (PHARE_WERROR_FLAGS -Wall -Wextra -pedantic -Werror)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
                        -Wno-unused-variable -Wno-unused-parameter
                        -Wno-gnu-zero-variadic-macro-arguments)
else()
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
                          -Wno-unused-variable -Wno-unused-parameter
                          -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)
endif()
endif()
