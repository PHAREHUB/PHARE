cmake_minimum_required (VERSION 3.20.1)
project(test-mhdsolver)
set(SOURCES test_mhd_solver.cpp)
add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  ${GTEST_INCLUDE_DIRS}
  )
target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  phare_amr
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_geom
  SAMRAI_hier
  SAMRAI_math
  SAMRAI_mesh
  SAMRAI_pdat
  SAMRAI_solv
  SAMRAI_tbox
  SAMRAI_xfer
  MPI::MPI_C
  ${GTEST_LIBS}
)
add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
