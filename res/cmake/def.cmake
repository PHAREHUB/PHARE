

if(testMPI)
  if (NOT DEFINED PHARE_MPI_PROCS)
    set(PHARE_MPI_PROCS 2)
  endif()
else()
  if (NOT DEFINED PHARE_MPI_PROCS)
    set(PHARE_MPI_PROCS 1)
  endif()
endif()

set (PHARE_HAS_HIGHFIVE "0")
if (HighFive)
 set (PHARE_HAS_HIGHFIVE "1")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE}")

# Pybind errors with clang, it is default in GCC
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsized-deallocation")
endif()

if(coverage AND NOT MSVC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -DHAVE_EXECINFO_H -g3")
endif()

set (PHARE_FLAGS ${PHARE_FLAGS})
set (PHARE_WERROR_FLAGS ${PHARE_FLAGS} ${PHARE_WERROR_FLAGS})

if(devMode)
if(MSVC)
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS} /W4 /WX)
else()
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS} -Wall -Wextra -pedantic -Werror
        -Wno-unused-variable -Wno-unused-parameter)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
      -Wno-gnu-zero-variadic-macro-arguments)
else() # !Clang
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
      -Wno-class-memaccess -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)
endif() # clang
endif() # msvc
endif() # devMode

function(phare_sanitize_ san cflags )
  set(CMAKE_REQUIRED_FLAGS ${san})
  check_cxx_compiler_flag( ${san} ADDRESS_SANITIZER)
  if (${ADDRESS_SANITIZER})
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${san} ${cflags}" PARENT_SCOPE)
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${san}"  PARENT_SCOPE)
    set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${san}"  PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Your compiler: ${CMAKE_CXX_COMPILER_ID} seems to not support ${san}")
  endif()
  unset(CMAKE_REQUIRED_FLAGS)
endfunction(phare_sanitize_)

if (asan)   # -Dasan=ON
  phare_sanitize_("-fsanitize=address" "-fno-omit-frame-pointer" )
endif()

if (ubsan)  # -Dubsan=ON
  phare_sanitize_("-fsanitize=undefined" "" )
endif()

# msan is not supported - it's not practical to configure - use valgrind

function(set_exe_paths_ binary)
  set_property(TEST ${binary}        PROPERTY ENVIRONMENT "PYTHONPATH=${PHARE_PYTHONPATH}")
  set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${LD_LIBRARY_PATH})
endfunction(set_exe_paths_)

function(add_phare_test_ binary directory)
  target_compile_options(${binary} PRIVATE ${PHARE_WERROR_FLAGS})
  set_exe_paths_(${binary})
  set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
  set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT PHARE_MPI_PROCS=${PHARE_MPI_PROCS})
endfunction(add_phare_test_)

function(add_no_mpi_phare_test binary directory)
  if(NOT testMPI OR (testMPI AND forceSerialTests))
    add_test(NAME ${binary} COMMAND ./${binary} WORKING_DIRECTORY ${directory})
    add_phare_test_(${binary} ${directory})
  endif()
endfunction(add_no_mpi_phare_test)

function(add_no_mpi_python3_test name file directory)
  if(NOT testMPI OR (testMPI AND forceSerialTests))
    add_test(NAME py3_${name} COMMAND python3 ${file} WORKING_DIRECTORY ${directory})
    set_exe_paths_(py3_${name})
  endif()
endfunction(add_no_mpi_python3_test)

if(testMPI)
  function(add_phare_test binary directory)
    add_test(NAME ${binary} COMMAND mpirun -n ${PHARE_MPI_PROCS} ./${binary} WORKING_DIRECTORY ${directory})
    add_phare_test_(${binary} ${directory})
  endfunction(add_phare_test)

  function(add_python3_test name file directory)
    add_test(NAME py3_${name} COMMAND mpirun -n ${PHARE_MPI_PROCS} python3 ${file} WORKING_DIRECTORY ${directory})
    set_exe_paths_(py3_${name})
  endfunction(add_python3_test)
else()
  function(add_phare_test binary directory)
    add_no_mpi_phare_test(${binary} ${directory})
  endfunction(add_phare_test)

  function(add_python3_test name file directory)
    add_no_mpi_python3_test(${name} ${file} ${directory})
  endfunction(add_python3_test)
endif(testMPI)
