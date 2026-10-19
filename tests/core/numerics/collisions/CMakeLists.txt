cmake_minimum_required (VERSION 3.20.1)

project(test-collisions)

foreach(target test-collisions test-collisions-relaxation)
  string(REPLACE "-" "_" source ${target})
  add_executable(${target} ${source}.cpp)

  target_include_directories(${target} PRIVATE
    ${GTEST_INCLUDE_DIRS}
    )

  target_link_libraries(${target} PRIVATE
    phare_initializer
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${target} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# relaxation runs, about a minute in debug, are kept out of the PR tier
if(TEST test-collisions-relaxation)
  set_property(TEST test-collisions-relaxation APPEND PROPERTY LABELS nightly)
endif()
