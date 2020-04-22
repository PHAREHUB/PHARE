cmake_minimum_required (VERSION 3.3)

project(test-GridLayout)

set(SOURCES_INC
  gridlayout_allocsize.h
  gridlayout_base_params.h
  gridlayout_cell_centered_coord.h
  gridlayout_deriv.h
  gridlayout_laplacian.h
  gridlayout_field_centered_coord.h
  gridlayout_indexing.h
  gridlayout_params.h
  gridlayout_test.h
  gridlayout_utilities.h
  test_linear_combinaisons_yee.h
   )

set(SOURCES_CPP
  gridlayout_amr.cpp
  gridlayout_allocsize.cpp
  gridlayout_cell_centered_coord.cpp
  gridlayout_deriv.cpp
  gridlayout_laplacian.cpp
  gridlayout_field_centered_coord.cpp
  gridlayout_indexing.cpp
  test_linear_combinaisons_yee.cpp
  test_nextprev.cpp
  test_main.cpp
   )
add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})



target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND "PYTHONPATH=${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/pyphare" ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gridlayout_generate_test_files.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/grid_yee.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

add_python3_test(grid_yee grid_yee.py ${CMAKE_CURRENT_BINARY_DIR})
