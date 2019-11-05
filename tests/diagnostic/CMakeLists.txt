cmake_minimum_required (VERSION 3.3)

project(test-diagnostics)
configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)

set(SOURCES_INC
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
  ${PROJECT_SOURCE_DIR}/include.h
  ${PROJECT_SOURCE_DIR}/func.h
  ${PROJECT_SOURCE_DIR}/defs.h
  ${PROJECT_SOURCE_DIR}/integrator.h
  ${PROJECT_SOURCE_DIR}/samrai_lifecycle.h
  ${PROJECT_SOURCE_DIR}/tag_strat.h
  ${PROJECT_SOURCE_DIR}/hierarchy.h
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
)

set(SOURCES_CPP
  ${PROJECT_SOURCE_DIR}/test_main.cpp
)



if(HighFive)
  add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

  target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
    )

  target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
  target_link_directories(${PROJECT_NAME} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )
  target_link_libraries(${PROJECT_NAME} PRIVATE
    phare_amr
    phare_solver
    phare_diagnostic
    gtest
    gmock
    ${HDF5_LIBRARIES})

  add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}
           WORKING_DIRECTORY ${PHARE_PROJECT_DIR})
  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})

  include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
endif()
