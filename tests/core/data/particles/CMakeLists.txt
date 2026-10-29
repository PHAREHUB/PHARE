cmake_minimum_required (VERSION 3.20.1)

project(test-particles)

function(_add_core_particles_test src_name test_name)

  add_executable(${test_name} ${src_name}.cpp)

  target_include_directories(${test_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${test_name} PRIVATE
    phare_core
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_core_particles_test)

_add_core_particles_test(test_particles test-particles)
_add_core_particles_test(test_particles_selecting test_particles_selecting)
_add_core_particles_test(test_particles_partitioning test_particles_partitioning)
_add_core_particles_test(test_particles_construction test_particles_construction)
_add_core_particles_test(test_particle_array_consistency test-particles-consistency)
_add_core_particles_test(test_particle_array_synchronisation test_particle_array_synchronisation)
_add_core_particles_test(test_particle_sorting test_particle_sorting)
_add_core_particles_test(test_interop test-particles-interop)
