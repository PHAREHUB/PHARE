cmake_minimum_required (VERSION 3.20.1)

project(test-ohm)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_initializer
  ${GTEST_LIBS})

add_custom_command(TARGET test-ohm
                   PRE_BUILD
                   COMMAND "PYTHONPATH=${PHARE_PYTHONPATH}" ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_ohm.py ${CMAKE_CURRENT_BINARY_DIR})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test_tileset_ohm test_tileset_ohm.cpp)

target_include_directories(test_tileset_ohm PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_directories(test_tileset_ohm PRIVATE
  ${HDF5_LIBRARY_PATH}
  )

target_link_libraries(test_tileset_ohm PRIVATE
  phare_core
  phare_simulator
  pybind11::embed
  ${GTEST_LIBS})

add_no_mpi_phare_test(test_tileset_ohm ${CMAKE_CURRENT_BINARY_DIR})
