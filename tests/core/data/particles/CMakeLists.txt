cmake_minimum_required (VERSION 3.9)

project(test-particles)

function(_particles_test src test_name)

  add_executable(${test_name} ${src})

  target_include_directories(${test_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${test_name} PRIVATE
    phare_core
    ${GTEST_LIBS})

  add_phare_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_particles_test)

_particles_test(test_main.cpp test-particles)
_particles_test(test_interop.cpp test-particles-interop)
