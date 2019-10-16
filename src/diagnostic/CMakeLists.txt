cmake_minimum_required (VERSION 3.3)

project(phare_diagnostic)

set(SOURCES_INC
  ${PROJECT_SOURCE_DIR}/diagnostic_manager.h
)

set(SOURCES_CPP
  ${PROJECT_SOURCE_DIR}/diagnostic_manager.cpp
)

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})

target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>/src/amr
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>/src/solver
  $<INSTALL_INTERFACE:include/phare/diagnostic>)

target_link_libraries(${PROJECT_NAME} PRIVATE initializer phare_solver)

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
