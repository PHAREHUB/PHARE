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

