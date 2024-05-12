

# Per compiler CXXFLAGS
set (PHARE_FLAGS ${PHARE_FLAGS} )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (PHARE_FLAGS ${PHARE_FLAGS} )
else() # !Clang
  set (PHARE_FLAGS ${PHARE_FLAGS} --param=min-pagesize=0 )
endif() # clang

set (PHARE_LINK_FLAGS )
set (PHARE_BASE_LIBS )

if(PGO_GEN)
  if(PGO_USE)
    message(FATAL_ERROR "cannot generate and use pgo at the same time.")
  endif()
  set (PHARE_FLAGS ${PHARE_FLAGS} -fprofile-generate -fprofile-update=prefer-atomic )
  set (PHARE_LINK_FLAGS "${PHARE_LINK_FLAGS} -fprofile-generate -fprofile-update=prefer-atomic" )
endif()

if(PGO_USE)
  set (PHARE_LINK_FLAGS ${PHARE_LINK_FLAGS} -fprofile-use )
  set (PHARE_FLAGS ${PHARE_FLAGS} -fprofile-use )
endif()


set (PHARE_WERROR_FLAGS ${PHARE_FLAGS} ${PHARE_WERROR_FLAGS})
set (PHARE_PYTHONPATH "${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/pyphare")
set (PHARE_MPIRUN_POSTFIX ${PHARE_MPIRUN_POSTFIX})

# now we see if we are running with configurator
if (phare_configurator)
  execute_process(
    COMMAND ./tools/config/cmake.sh "${CMAKE_COMMAND}" "${CMAKE_CXX_COMPILER}" "${Python_EXECUTABLE}"
    WORKING_DIRECTORY ${PHARE_PROJECT_DIR}
    COMMAND_ERROR_IS_FATAL ANY
  )
  include("${PHARE_PROJECT_DIR}/tools/config/local.cmake")
endif(phare_configurator)

# Link Time Optimisation flags - is disabled if coverage is enabled
set (PHARE_INTERPROCEDURAL_OPTIMIZATION FALSE)
if(withIPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT PHARE_INTERPROCEDURAL_OPTIMIZATION OUTPUT error)
endif(withIPO)


set (PHARE_WITH_CCACHE FALSE)
if(devMode) # -DdevMode=ON

  # Having quotes on strings here has lead to quotes being added to the compile string, so avoid.

  set (_Werr ${PHARE_WERROR_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-unused-variable -Wno-unused-parameter)
  set (_Werr ${_Werr} -Wdouble-promotion -Wuninitialized )

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (_Werr ${_Werr} -Wno-gnu-zero-variadic-macro-arguments)

  else() # !Clang
    set (_Werr ${_Werr} -Wno-class-memaccess -Wno-unused-but-set-variable -Wno-unused-but-set-parameter)

  endif() # clang

  set (PHARE_WERROR_FLAGS ${_Werr})

  if(withCcache)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
      set(PHARE_WITH_CCACHE TRUE)
    endif()
  endif()
endif(devMode)

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
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    phare_sanitize_("-fsanitize=address" "-fno-omit-frame-pointer" )
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    phare_sanitize_("-fsanitize=address -shared-libsan" "-fno-omit-frame-pointer" )
  else()
    message(FATAL_ERROR "ASAN Unhandled compiler: ${CMAKE_CXX_COMPILER_ID}")
  endif()
  set(testDuringBuild OFF) # can need LD_PRELOAD/etc
endif(asan)

if (ubsan)  # -Dubsan=ON
  phare_sanitize_("-fsanitize=undefined" "" )
endif(ubsan)

# msan is not supported - it's not practical to configure - use valgrind

# test functions below

# public test functions
#
#   add_phare_test($binary $directory)
#    execute binary in target directory, with mpirun when -DtestMPI=ON
#
#   add_python3_test($name $file $directory)
#    launch python3 file described by name in target directory, with mpirun when -DtestMPI=ON
#
#   add_no_mpi_phare_test($binary $directory)
#    execute binary in target directory, does not run when -DtestMPI=ON
#
#   add_no_mpi_python3_test($name $file $directory)
#    launch python3 file described by name in target directory, does not run when -DtestMPI=ON
#
#  phare_exec(level target exe directory)
#    execute exe identified by target in directory
#      if level >= PHARE_EXEC_LEVEL_MIN AND level <= PHARE_EXEC_LEVEL_MAX
#
#  phare_python3_exec(level target file directory $ARGS)
#    execute file identified by target in directory
#      if level >= PHARE_EXEC_LEVEL_MIN AND level <= PHARE_EXEC_LEVEL_MAX
#
#  Note to developers - do not use cmake variable function arguments for functions
#      phare_python3_exec
#      phare_mpi_python3_exec
#   if these function calls are to files executing python unit tests as they will interfere

if (test AND ${PHARE_EXEC_LEVEL_MIN} GREATER 0) # 0 = no tests

  if (NOT DEFINED PHARE_MPI_PROCS)
    set(PHARE_MPI_PROCS 1)
    if(testMPI)
      set(PHARE_MPI_PROCS 2)
    endif()
  endif()

  function(set_exe_paths_ binary)
    set_property(TEST ${binary}        PROPERTY ENVIRONMENT "PYTHONPATH=${PHARE_PYTHONPATH}")
    # ASAN detects leaks by default, even in system/third party libraries
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT PHARE_SKIP_CLI=1 )
  endfunction(set_exe_paths_)

  function(add_phare_test_ binary directory)
    target_compile_options(${binary} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
    set_exe_paths_(${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT PHARE_MPI_PROCS=${PHARE_MPI_PROCS})
    if(testDuringBuild)
      add_custom_command(
           TARGET ${binary}
           POST_BUILD
           WORKING_DIRECTORY ${directory}
           COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> -R \"^${binary}$$\" --output-on-failures
      )
    endif(testDuringBuild)
  endfunction(add_phare_test_)

  function(add_no_mpi_phare_test binary directory)
    if(NOT testMPI OR (testMPI AND forceSerialTests))
      add_test(NAME ${binary} COMMAND ./${binary} WORKING_DIRECTORY ${directory})
      add_phare_test_(${binary} ${directory})
    else()
      # this prevents building targets even when added via "add_executable"
      set_target_properties(${binary} PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
    endif()
  endfunction(add_no_mpi_phare_test)

  function(add_no_mpi_python3_test name file directory)
    if(NOT testMPI OR (testMPI AND forceSerialTests))
      add_test(NAME py3_${name} COMMAND python3 -u ${file} WORKING_DIRECTORY ${directory})
      set_exe_paths_(py3_${name})
    endif()
  endfunction(add_no_mpi_python3_test)

  if(testMPI)
    function(add_phare_test binary directory)
      add_test(NAME ${binary} COMMAND mpirun -n ${PHARE_MPI_PROCS} ${PHARE_MPIRUN_POSTFIX} ./${binary} WORKING_DIRECTORY ${directory})
      add_phare_test_(${binary} ${directory})
    endfunction(add_phare_test)

    function(add_python3_test name file directory)
      add_test(NAME py3_${name} COMMAND mpirun -n ${PHARE_MPI_PROCS} ${PHARE_MPIRUN_POSTFIX} python3 -u ${file} WORKING_DIRECTORY ${directory})
      set_exe_paths_(py3_${name})
    endfunction(add_python3_test)

    function(add_mpi_python3_test N name file directory)
      add_test(NAME py3_${name}_mpi_n_${N} COMMAND mpirun -n ${N} ${PHARE_MPIRUN_POSTFIX} python3 ${file} WORKING_DIRECTORY ${directory})
      set_exe_paths_(py3_${name}_mpi_n_${N})
    endfunction(add_mpi_python3_test)

  else()
    function(add_phare_test binary directory)
      add_no_mpi_phare_test(${binary} ${directory})
    endfunction(add_phare_test)

    function(add_python3_test name file directory)
      add_no_mpi_python3_test(${name} ${file} ${directory})
    endfunction(add_python3_test)

    function(add_mpi_python3_test N name file directory)
      # do nothing
    endfunction(add_mpi_python3_test)
  endif(testMPI)


  function(add_phare_build_test binary directory)
  endfunction(add_phare_build_test)


  if(DEFINED GTEST_ROOT)
    set(GTEST_ROOT ${GTEST_ROOT} CACHE PATH "Path to googletest")
    find_package(GTest REQUIRED)
    set(GTEST_LIBS GTest::GTest GTest::Main)
  else()
    set(GTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

    if (NOT EXISTS ${GTEST_ROOT})
      execute_process(COMMAND ${Git} clone https://github.com/google/googletest ${GTEST_ROOT})
    endif()

    add_subdirectory(subprojects/googletest)
    set(GTEST_INCLUDE_DIRS
      $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>)
    set(GTEST_LIBS gtest gmock)

  endif()

  function(phare_exec level target exe directory)
    if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      add_test(NAME ${target} COMMAND ${exe} WORKING_DIRECTORY ${directory})
    endif()
  endfunction(phare_exec)
  # use
  #  phare_exec(1 test_id ./binary ${CMAKE_CURRENT_BINARY_DIR})


  function(phare_mpi_python3_exec level N target file directory)
    if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      string (REPLACE ";" " " CLI_ARGS "${ARGN}")
      if(${N} EQUAL 1)
        add_test(
            NAME py3_${target}
            COMMAND python3 -u ${file} ${CLI_ARGS}
            WORKING_DIRECTORY ${directory})
        set_exe_paths_(py3_${target})
      else()
        add_test(
            NAME py3_${target}_mpi_n_${N}
            COMMAND mpirun -n ${N} ${PHARE_MPIRUN_POSTFIX} python3 -u ${file} ${CLI_ARGS}
            WORKING_DIRECTORY ${directory})
        set_exe_paths_(py3_${target}_mpi_n_${N})
      endif()
    endif()
  endfunction(phare_mpi_python3_exec)
  # use
  #  phare_mpi_python3_exec(1 2 test_id script.py ${CMAKE_CURRENT_BINARY_DIR} $ARGS)


  if(testMPI)
    function(phare_python3_exec level target file directory)
      phare_mpi_python3_exec(${level} 2 ${target} ${file} ${directory})
    endfunction(phare_python3_exec)
  else()
    function(phare_python3_exec level target file directory)
      if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
        string (REPLACE ";" " " CLI_ARGS "${ARGN}")
        add_test(NAME py3_${target} COMMAND python3 -u ${file} ${CLI_ARGS} WORKING_DIRECTORY ${directory})
        set_exe_paths_(py3_${target})
      endif()
    endfunction(phare_python3_exec)
  endif(testMPI)
  # use
  #  phare_python3_exec(1 test_id script.py ${CMAKE_CURRENT_BINARY_DIR} $ARGS)


  set(GTEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIRS} ${PHARE_PROJECT_DIR})

  enable_testing()

endif()

# useful to see what's available after importing a package
function(phare_print_all_vars)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endfunction(phare_print_all_vars)

