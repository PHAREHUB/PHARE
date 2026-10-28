cmake_minimum_required (VERSION 3.20.1)

project(test-models)



set(SOURCES_CPP
  test_models.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  phare_initializer
  ${GTEST_LIBS})


add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



add_executable(test-external-field-manager test_external_field_manager.cpp)

target_include_directories(test-external-field-manager PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(test-external-field-manager PRIVATE
  phare_amr
  phare_initializer
  ${GTEST_LIBS})

add_no_mpi_phare_test(test-external-field-manager ${CMAKE_CURRENT_BINARY_DIR})
