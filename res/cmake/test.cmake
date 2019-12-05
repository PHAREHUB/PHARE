
if (test)

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
  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/field)
  add_subdirectory(tests/initializer)
  add_subdirectory(tests/diagnostic)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ions)
  add_subdirectory(tests/core/data/ion_population)
  add_subdirectory(tests/core/data/maxwellian_particle_initializer)
  add_subdirectory(tests/core/data/particle_initializer)
  add_subdirectory(tests/core/utilities/box)
  add_subdirectory(tests/core/utilities/particle_selector)
  add_subdirectory(tests/core/utilities/partitionner)
  add_subdirectory(tests/core/utilities/range)
  add_subdirectory(tests/core/utilities/index)
  add_subdirectory(tests/core/numerics/boundary_condition)
  add_subdirectory(tests/core/numerics/interpolator)
  add_subdirectory(tests/core/numerics/pusher)
  add_subdirectory(tests/core/numerics/ampere)
  add_subdirectory(tests/core/numerics/faraday)
  add_subdirectory(tests/core/numerics/ohm)

endif()
