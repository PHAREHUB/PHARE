



if (test AND ${PHARE_EXEC_LEVEL_MIN} GREATER 0) # 0 = no tests

  configure_file(${CMAKE_SOURCE_DIR}/tests/__init__.py ${CMAKE_BINARY_DIR}/tests/__init__.py @ONLY)


  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/grid)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ions)
  add_subdirectory(tests/core/data/electrons)
  add_subdirectory(tests/core/data/ion_population)
  add_subdirectory(tests/core/data/maxwellian_particle_initializer)
  add_subdirectory(tests/core/data/particle_initializer)
  add_subdirectory(tests/core/data/mhd_state)
  add_subdirectory(tests/core/utilities/box)
  add_subdirectory(tests/core/utilities/range)
  add_subdirectory(tests/core/utilities/index)
  add_subdirectory(tests/core/utilities/indexer)
  add_subdirectory(tests/core/utilities/cellmap)
  #add_subdirectory(tests/core/numerics/boundary_condition)
  add_subdirectory(tests/core/numerics/interpolator)
  add_subdirectory(tests/core/numerics/pusher)
  add_subdirectory(tests/core/numerics/ampere)
  add_subdirectory(tests/core/numerics/faraday)
  add_subdirectory(tests/core/numerics/ohm)
  add_subdirectory(tests/core/numerics/ion_updater)
  add_subdirectory(tests/core/numerics/mock_mhd_simulator)


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
  add_subdirectory(tests/functional/dispersion)
  add_subdirectory(tests/functional/ionIonBeam)
  add_subdirectory(tests/functional/conservation)
  add_subdirectory(tests/functional/harris)

  add_subdirectory(pyphare/pyphare_tests/test_pharesee/)
  add_subdirectory(pyphare/pyphare_tests/pharein/)
  add_subdirectory(pyphare/pyphare_tests/test_core/)


endif()
