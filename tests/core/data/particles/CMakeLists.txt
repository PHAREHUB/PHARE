cmake_minimum_required (VERSION 3.20.1)

project(test_particles)

function(_add_core_particles_test src_name)

  add_executable(${src_name} ${src_name}.cpp)

  target_include_directories(${src_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${src_name} PRIVATE
    phare_initializer
    phare_core
    ${GTEST_LIBS}
  )

  add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_core_particles_test)

_add_core_particles_test(test_particles)
_add_core_particles_test(test_particles_selecting)
_add_core_particles_test(test_particle_array_consistency)
