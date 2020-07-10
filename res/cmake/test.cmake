
# public functions
#
#   add_phare_test($project $directory)
#
#
#   add_python3_test($name $file $directory)
#
#
#   add_no_mpi_phare_test($project $directory)
#
#
#   add_no_mpi_python3_test($name $file $directory)
#
#

if (test)

  function(set_exe_paths_ binary)
    set_property(TEST ${binary}        PROPERTY ENVIRONMENT "PYTHONPATH=${PHARE_PYTHONPATH}")
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")
  endfunction(set_exe_paths_)

  function(add_phare_test_ binary directory)
    target_compile_options(${binary} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
    set_exe_paths_(${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT GMON_OUT_PREFIX=gprof.${binary})
    set_property(TEST ${binary} APPEND PROPERTY ENVIRONMENT PHARE_MPI_PROCS=${PHARE_MPI_PROCS})
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


  if(DEFINED GTEST_ROOT)
    set(GTEST_ROOT ${GTEST_ROOT} CACHE PATH "Path to googletest")
    find_package(GTest REQUIRED)
    set(GTEST_LIBS GTest::GTest GTest::Main)
  else()
    set(GTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

    if (NOT EXISTS ${GTEST_ROOT})
      execute_process(
        COMMAND ${Git} clone https://github.com/google/googletest ${GTEST_ROOT}
      )
    endif()

    add_subdirectory(subprojects/googletest)
    set(GTEST_INCLUDE_DIRS
      $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>)
    set(GTEST_LIBS gtest gmock)

  endif()

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
  add_subdirectory(tests/amr/periodicity)
  add_subdirectory(tests/amr/models)
  add_subdirectory(tests/amr/multiphysics_integrator)


  add_subdirectory(tests/diagnostic)


  add_subdirectory(tests/simulator)


  add_subdirectory(pyphare/pyphare_tests/test_pharesee/)
  add_subdirectory(pyphare/pyphare_tests/pharein/)
  add_subdirectory(pyphare/pyphare_tests/test_core/)


endif()

