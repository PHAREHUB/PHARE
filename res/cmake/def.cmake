

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

if(MSVC)
  set (PHARE_WERROR_FLAGS /W4 /WX)
else()
  set (PHARE_WERROR_FLAGS -Wall -Wextra -pedantic -Werror
        -Wno-unused-variable -Wno-unused-parameter)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
      -Wno-gnu-zero-variadic-macro-arguments)
else() # !Clang
  set (PHARE_WERROR_FLAGS ${PHARE_WERROR_FLAGS}
      -Wno-class-memaccess -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)
endif() # clang
endif() # msvc

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

function(add_phare_test_ binary directory)
  target_compile_options(${binary} PRIVATE ${PHARE_WERROR_FLAGS})
  set_tests_properties(${binary} PROPERTIES ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
  set_tests_properties(${binary} PROPERTIES ENVIRONMENT
      PYTHONPATH=${CMAKE_BINARY_DIR})
endfunction(add_phare_test_)

if(testMPI)
  function(add_phare_test binary directory)
    add_test(NAME ${binary} COMMAND mpirun -n ${PHARE_MPI_PROCS} ./${binary} WORKING_DIRECTORY ${directory})
    add_phare_test_(${binary} ${directory})
  endfunction(add_phare_test)
else()
  function(add_phare_test binary directory)
    add_test(NAME ${binary} COMMAND ./${binary} WORKING_DIRECTORY ${directory})
    add_phare_test_(${binary} ${directory})
  endfunction(add_phare_test)
endif(testMPI)
