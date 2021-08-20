
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
#  phare_python3_exec(level target file directory)
#    execute file identified by target in directory
#      if level >= PHARE_EXEC_LEVEL_MIN AND level <= PHARE_EXEC_LEVEL_MAX
#

# this define breaks pgcc / and appears unnecessary
# sed -i 's/# define GTEST_CAN_STREAM_RESULTS_ 1/# define GTEST_CAN_STREAM_RESULTS_ 0/g' subprojects/googletest/googletest/include/gtest/internal/gtest-port.h

if (test AND ${PHARE_EXEC_LEVEL_MIN} GREATER 0) # 0 = no tests

  if (NOT DEFINED PHARE_MPI_PROCS)
    set(PHARE_MPI_PROCS 1)
    if(testMPI)
      set(PHARE_MPI_PROCS 2)
    endif()
  endif()

  function(set_exe_paths_ binary)
    set_property(TEST ${binary}        PROPERTY ENVIRONMENT "PYTHONPATH=${PHARE_PYTHONPATH}")
  endfunction(set_exe_paths_)

  function(add_phare_test_ binary)
    target_compile_options(${binary} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
    target_link_options(${binary} PRIVATE ${PHARE_BIN_FLAGS})
    set_exe_paths_(${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT PHARE_MPI_PROCS=${PHARE_MPI_PROCS})
  endfunction(add_phare_test_)

  function(add_no_mpi_phare_test binary directory)
    if(NOT testMPI OR (testMPI AND forceSerialTests))
      add_test(NAME ${binary} COMMAND ./${binary} WORKING_DIRECTORY ${directory})
      add_phare_test_(${binary})
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
      add_test(NAME ${binary} COMMAND mpirun -n ${PHARE_MPI_PROCS} ./${binary} WORKING_DIRECTORY ${directory})
      add_phare_test_(${binary})
    endfunction(add_phare_test)

    function(add_python3_test name file directory)
      add_test(NAME py3_${name} COMMAND mpirun -n ${PHARE_MPI_PROCS} python3 -u ${file} WORKING_DIRECTORY ${directory})
      set_exe_paths_(py3_${name})
    endfunction(add_python3_test)

    function(add_mpi_python3_test N name file directory)
      add_test(NAME py3_${name}_mpi_n_${N} COMMAND mpirun -n ${N} python3 ${file} WORKING_DIRECTORY ${directory})
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


  function(phare_exec level target exe directory)
    if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      add_test(NAME ${target} COMMAND ${exe} WORKING_DIRECTORY ${directory})
    endif()
  endfunction(phare_exec)
  # use
  #  phare_exec(1 test_id ./binary ${CMAKE_CURRENT_BINARY_DIR})

  function(phare_python3_exec level target file directory)
    if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      add_test(NAME py3_${target} COMMAND python3 -u ${file} WORKING_DIRECTORY ${directory})
      set_exe_paths_(py3_${target})
    endif()
  endfunction(phare_python3_exec)
  # use
  #  phare_python3_exec(1 test_id script.py ${CMAKE_CURRENT_BINARY_DIR})


  function(phare_mpi_python3_exec level N target file directory)
    if(${level} GREATER_EQUAL ${PHARE_EXEC_LEVEL_MIN} AND ${level} LESS_EQUAL ${PHARE_EXEC_LEVEL_MAX})
      if(${N} EQUAL 1)
        add_test(NAME py3_${target} COMMAND python3 -u ${file} WORKING_DIRECTORY ${directory})
        set_exe_paths_(py3_${target})
      else()
        add_test(NAME py3_${target}_mpi_n_${N} COMMAND mpirun -n ${N} python3 -u ${file} WORKING_DIRECTORY ${directory})
        set_exe_paths_(py3_${target}_mpi_n_${N})
      endif()
    endif()
  endfunction(phare_mpi_python3_exec)
  # use
  #  phare_mpi_python3_exec(1 2 test_id script.py ${CMAKE_CURRENT_BINARY_DIR})

  set(GTEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIRS} ${PHARE_PROJECT_DIR})

  enable_testing()

  configure_file(${CMAKE_SOURCE_DIR}/tests/__init__.py ${CMAKE_BINARY_DIR}/tests/__init__.py @ONLY)


  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/field)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ions)
  add_subdirectory(tests/core/data/electrons)
  add_subdirectory(tests/core/data/ion_population)
  add_subdirectory(tests/core/data/maxwellian_particle_initializer)
  add_subdirectory(tests/core/data/particle_initializer)
  add_subdirectory(tests/core/utilities/box)
  add_subdirectory(tests/core/utilities/partitionner)
  add_subdirectory(tests/core/utilities/range)
  add_subdirectory(tests/core/utilities/index)
  add_subdirectory(tests/core/numerics/boundary_condition)
  add_subdirectory(tests/core/numerics/interpolator)
  add_subdirectory(tests/core/numerics/pusher)
  add_subdirectory(tests/core/numerics/ampere)
  add_subdirectory(tests/core/numerics/faraday)
  add_subdirectory(tests/core/numerics/ohm)
  add_subdirectory(tests/core/numerics/ion_updater)


  add_subdirectory(tests/initializer)


  add_subdirectory(tests/amr/data/particles)
  add_subdirectory(tests/amr/data/field/coarsening)
  add_subdirectory(tests/amr/data/field/copy_pack)
  add_subdirectory(tests/amr/data/field/geometry)
  add_subdirectory(tests/amr/data/field/overlap)
  add_subdirectory(tests/amr/data/field/refine)
  add_subdirectory(tests/amr/data/field/variable)
  add_subdirectory(tests/amr/data/field/time_interpolate)
  add_subdirectory(tests/amr/resources_manager)
  add_subdirectory(tests/amr/messengers)
  add_subdirectory(tests/amr/models)
  add_subdirectory(tests/amr/multiphysics_integrator)
  add_subdirectory(tests/amr/tagging)

  add_subdirectory(tests/diagnostic)


  add_subdirectory(tests/simulator)

  add_subdirectory(tests/functional/alfven_wave)
  add_subdirectory(tests/functional/td)
  add_subdirectory(tests/functional/translation)
  add_subdirectory(tests/functional/tdtagged)
  add_subdirectory(tests/functional/shock)
  add_subdirectory(tests/functional/conservation)
  add_subdirectory(tests/functional/harris)

  add_subdirectory(pyphare/pyphare_tests/test_pharesee/)
  add_subdirectory(pyphare/pyphare_tests/pharein/)
  add_subdirectory(pyphare/pyphare_tests/test_core/)


endif()
