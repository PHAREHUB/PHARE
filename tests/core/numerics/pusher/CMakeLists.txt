cmake_minimum_required (VERSION 3.20.1)

project(test-pusher)

function(_add_pusher_test name src)

  add_executable(${name} ${src})

  target_include_directories(${name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
    )

  target_link_libraries(${name} PRIVATE
    phare_core
    phare_initializer
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_pusher_test)

_add_pusher_test(${PROJECT_NAME} test_pusher.cpp)
_add_pusher_test(test-particle-precision test_particle_precision.cpp)

add_custom_command(TARGET test-pusher
                   PRE_BUILD
                   COMMAND "PYTHONPATH=${PHARE_PYTHONPATH}" ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_pusher.py ${CMAKE_CURRENT_BINARY_DIR})
