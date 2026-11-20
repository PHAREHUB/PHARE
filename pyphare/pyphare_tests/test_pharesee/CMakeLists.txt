cmake_minimum_required (VERSION 3.20.1)

project(test-pharesee)


add_python3_test(test-pharesee-geometry_1d test_geometry.py ${PROJECT_SOURCE_DIR})
add_python3_test(test-pharesee-geometry_2d test_geometry_2d.py ${PROJECT_SOURCE_DIR})

#phare_mpi_python3_exec(11 2 test-pharesee-hierarchy test_hierarchy.py ${PROJECT_SOURCE_DIR})

if(NOT lowResourceTests)

# gcov takes 30 hours
add_python3_test(test-pharesee-hierarchy test_hierarchy.py ${PROJECT_SOURCE_DIR})

endif()

# Sample runs of the pharesee.amr tests (amr/gen_sample_data.py): small simulations written
# in phareh5 AND vtkhdf, under ${CMAKE_BINARY_DIR}/sample_data. One ctest fixture per
# (case, format), so "ctest -j" generates them in parallel, each with its own timeout (a run
# already up to date is skipped). A test file requires the fixtures of the cases it reads
# (PHARESEE_AMR_CASES_<file stem>). Those tests never skip: missing data is an error
# (amr/sample_data.py require()). So:
#  - lowResourceTests, coverage, asan: no generation (the simulations run for hours under gcov,
#    > 20 min under asan), the files reading the sample runs are not registered; these tests
#    are pure Python, so the other builds cover them;
#  - PHARE_EXEC_LEVEL_MIN > 1: add_python3_test registers nothing, so nothing is generated;
#  - the 3-rank run (${CMAKE_BINARY_DIR}/sample_data_mpi3) needs mpirun: generated, and
#    amr/test_mpi3.py (its only reader) registered, only with testMPI.
set(PHARESEE_SAMPLE_DATA ${CMAKE_BINARY_DIR}/sample_data)
set(PHARESEE_SAMPLE_CASES 1d 2d_static 2d_tagged 3d 2d_mhd)
set(PHARESEE_SAMPLE_FORMATS h5 vtk)
# sample cases read by each amr test file (both formats of a case: sample_data.require)
set(PHARESEE_AMR_CASES_test_acceptance_issues 2d_tagged)
set(PHARESEE_AMR_CASES_test_compat 1d 2d_static 2d_tagged 3d 2d_mhd)
set(PHARESEE_AMR_CASES_test_differential_h5 2d_static)
set(PHARESEE_AMR_CASES_test_field 2d_tagged)
set(PHARESEE_AMR_CASES_test_io_phareh5 1d 2d_static 2d_tagged 3d)
set(PHARESEE_AMR_CASES_test_io_vtkhdf 1d 2d_static 2d_tagged 3d)
set(PHARESEE_AMR_CASES_test_mhd 2d_mhd)
set(PHARESEE_AMR_CASES_test_plot 2d_tagged)
set(PHARESEE_AMR_CASES_test_run_amr 1d 2d_static 2d_tagged 3d)
set(PHARESEE_AMR_SAMPLE_TESTS
    test_acceptance_issues test_compat test_differential_h5 test_field test_io_phareh5
    test_io_vtkhdf test_mhd test_plot test_run_amr)
set(PHARESEE_AMR_SKIPPED_TESTS test_mpi3)  # files not registered by the loop below
if(lowResourceTests OR coverage OR asan)
  list(APPEND PHARESEE_AMR_SKIPPED_TESTS ${PHARESEE_AMR_SAMPLE_TESTS})
elseif(NOT ${PHARE_EXEC_LEVEL_MIN} GREATER ${PHARE_UNLEVELED_TEST_LEVEL})
  # timeout of each generation: >= 4x its measured time, at least 600 s. Measured locally:
  # 1d 3 s, 3d 3 s, 2d_static 12 s, 2d_tagged 136 s, 2d_mhd 2 s (per format; the CI agents
  # are ~2.5x slower)
  set(PHARESEE_TIMEOUT_1d 600)
  set(PHARESEE_TIMEOUT_2d_static 600)
  set(PHARESEE_TIMEOUT_2d_tagged 1200)
  set(PHARESEE_TIMEOUT_3d 600)
  set(PHARESEE_TIMEOUT_2d_mhd 600)
  foreach(case ${PHARESEE_SAMPLE_CASES})
    foreach(fmt ${PHARESEE_SAMPLE_FORMATS})
      set(gen_test py3_test-pharesee-amr-sample-data-${case}-${fmt})
      add_test(NAME ${gen_test}
               COMMAND ${Python_EXECUTABLE} -u amr/gen_sample_data.py --root ${PHARESEE_SAMPLE_DATA}
                       ${case} ${fmt}
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
      set_exe_paths_(${gen_test})
      set_tests_properties(${gen_test} PROPERTIES FIXTURES_SETUP pharesee_sample_data_${case}_${fmt}
                           TIMEOUT ${PHARESEE_TIMEOUT_${case}})
    endforeach()
  endforeach()
  if(testMPI)
    foreach(fmt h5 vtk)
      set(gen_test py3_test-pharesee-amr-sample-data-mpi3-${fmt})
      add_test(NAME ${gen_test}
               COMMAND mpirun -n 3 ${PHARE_MPIRUN_POSTFIX} python3 -u amr/gen_sample_data.py
                       --root ${PHARESEE_SAMPLE_DATA} --mpi3 ${fmt}
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
      set_exe_paths_(${gen_test})
      set_test_procs_(${gen_test} 3)
      set_tests_properties(${gen_test} PROPERTIES FIXTURES_SETUP pharesee_sample_data_mpi3
                           TIMEOUT 1800)
    endforeach()
    set(amr_ctest py3_test-pharesee-amr-test_mpi3)
    add_python3_test(test-pharesee-amr-test_mpi3 amr/test_mpi3.py ${PROJECT_SOURCE_DIR})
    if(TEST ${amr_ctest})
      set_property(TEST ${amr_ctest} APPEND PROPERTY ENVIRONMENT PHARESEE_SAMPLE_DATA=${PHARESEE_SAMPLE_DATA})
      set_tests_properties(${amr_ctest} PROPERTIES FIXTURES_REQUIRED pharesee_sample_data_mpi3)
    endif()
  endif()
endif()

# pharesee.amr (new reader): pure Python unittest files, one ctest each, named
# py3_test-pharesee-amr-<file stem>. CONFIGURE_DEPENDS: a new amr/test_*.py file is picked up
# at the next build.
file(GLOB PHARESEE_AMR_TESTS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/amr/test_*.py)
foreach(amr_test ${PHARESEE_AMR_TESTS})
  get_filename_component(amr_test_name ${amr_test} NAME_WE)
  if(amr_test_name IN_LIST PHARESEE_AMR_SKIPPED_TESTS)
    continue()
  endif()
  set(amr_ctest py3_test-pharesee-amr-${amr_test_name})
  add_python3_test(test-pharesee-amr-${amr_test_name} amr/${amr_test_name}.py ${PROJECT_SOURCE_DIR})
  if(TEST ${amr_ctest} AND amr_test_name IN_LIST PHARESEE_AMR_SAMPLE_TESTS)
    set_property(TEST ${amr_ctest} APPEND PROPERTY ENVIRONMENT PHARESEE_SAMPLE_DATA=${PHARESEE_SAMPLE_DATA})
    set(fixtures)
    foreach(case ${PHARESEE_AMR_CASES_${amr_test_name}})
      foreach(fmt ${PHARESEE_SAMPLE_FORMATS})
        list(APPEND fixtures pharesee_sample_data_${case}_${fmt})
      endforeach()
    endforeach()
    set_tests_properties(${amr_ctest} PROPERTIES FIXTURES_REQUIRED "${fixtures}")
  endif()
endforeach()
