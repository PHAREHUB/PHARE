cmake_minimum_required (VERSION 3.20.1)

project(test-simulator)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "utilities" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

phare_python3_exec(9 validation           test_validation.py      ${CMAKE_CURRENT_BINARY_DIR})
phare_python3_exec(9 data-wrangler        data_wrangler.py        ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(periodicity test_init_periodicity.py ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(mhd_only_api test_mhd_only_api.py ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(time_step test_time_step.py ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(patch_overlap test_patch_overlap.py ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(multi_species_momentum test_multi_species_momentum.py ${CMAKE_CURRENT_BINARY_DIR})

if (useExceptionsInsteadOfMPIAbort)
  add_no_mpi_python3_test(exceptions test_exceptions.py ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  phare_python3_exec(9       diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       restarts    test_restarts.py     ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       run         test_run.py          ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       vtk_diagnostics test_vtk_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       seed_tiling test_seed_tiling.py  ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       mhd_shared_faces test_mhd_shared_faces.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_python3_exec(9       fully_kinetic test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_python3_exec(9       fully_kinetic_2d test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1)
  phare_python3_exec(10      fully_kinetic_2d_o23 test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 2 3)
  phare_python3_exec(9       fully_kinetic_3d test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 3 1)
  phare_python3_exec(10      fully_kinetic_3d_o23 test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 3 2 3)
  # current filter: J ghosts at patch corners, Gauss with the filtered charge
  phare_python3_exec(9       fully_kinetic_2d_filter test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1)
  phare_python3_exec(10      fully_kinetic_3d_filter test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 3 1 --filter-passes 1)
  # Marder cleaning: the start-up Gauss residual decays as predicted, mode by mode
  phare_python3_exec(9       fully_kinetic_cleaning test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 1 --cleaning 0.5)
  phare_python3_exec(9       fully_kinetic_2d_cleaning test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --cleaning 0.5)
  phare_python3_exec(9       fully_kinetic_2d_filter_cleaning test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1 --cleaning 0.5)
  # one static refined box (P4 T1, T2): fine time step, finite E, root Gauss drift outside the box
  phare_python3_exec(9       fully_kinetic_amr test_fully_kinetic_amr.py ${CMAKE_CURRENT_BINARY_DIR} 1 1 --cleaning 0)
  phare_python3_exec(9       fully_kinetic_amr_2d test_fully_kinetic_amr.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1 --cleaning 0.5)
  # P5a, T15: tagging on moving light pulses, level 1 regridded every 2 coarse steps
  phare_python3_exec(9       fully_kinetic_amr_regrid test_fully_kinetic_amr_regrid.py ${CMAKE_CURRENT_BINARY_DIR} 1 1 --cleaning 0.5)
  phare_python3_exec(9       fully_kinetic_amr_regrid_2d test_fully_kinetic_amr_regrid.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1 --cleaning 0.5)
  # T9 residual: a seeded B_z in vacuum does not grow with a refined box (~40 s, Debug)
  phare_python3_exec(10      fully_kinetic_amr_vacuum_stability test_fully_kinetic_amr_vacuum_stability.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_test_labels(fully_kinetic_amr_vacuum_stability physics)
  # phare_python3_exec(9       tagging     test_tagging.py      ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2

  if(testMPI)
    phare_mpi_python3_exec(9 3 diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 4 diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR})

    # doesn't make sense in serial
    phare_mpi_python3_exec(9 3 load_balancing test_load_balancing.py ${CMAKE_CURRENT_BINARY_DIR})

    phare_mpi_python3_exec(9 4 vtk_diagnostics test_vtk_diagnostics.py ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 3 seed_tiling test_seed_tiling.py ${CMAKE_CURRENT_BINARY_DIR})
    # patch border sums and particle exchange of the fully kinetic model between ranks
    phare_mpi_python3_exec(9 3 fully_kinetic test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 3 fully_kinetic_2d test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1)
    phare_mpi_python3_exec(9 3 fully_kinetic_3d test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 3 1)
    phare_mpi_python3_exec(9 3 fully_kinetic_2d_filter test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1)
    phare_mpi_python3_exec(9 3 fully_kinetic_2d_filter_cleaning test_fully_kinetic.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1 --cleaning 0.5)
    phare_mpi_python3_exec(9 3 fully_kinetic_amr_2d test_fully_kinetic_amr.py ${CMAKE_CURRENT_BINARY_DIR} 2 1 --filter-passes 1 --cleaning 0.5)
    phare_mpi_python3_exec(9 4 mhd_shared_faces test_mhd_shared_faces.py ${CMAKE_CURRENT_BINARY_DIR})

  else() # serial only
    phare_python3_exec(11 test_diagnostic_timestamps test_diagnostic_timestamps.py ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)


endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)

add_subdirectory(initialize)
add_subdirectory(refinement)
add_subdirectory(advance)
add_subdirectory(amr_convergence)

