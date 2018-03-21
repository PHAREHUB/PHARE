cmake_minimum_required(VERSION 3.3)
project(phare_amr)

set( SOURCES_INC
     tools/resources_manager.h
     tools/resources_manager_utilities.h
     tools/resources_guards.h
   )
set( SOURCES_CPP
     tools/resources_manager.cpp
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/amr>)

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${SAMRAI_INCLUDE_DIRS}>)


target_link_libraries(${PROJECT_NAME} PUBLIC phare_core)
target_link_libraries(${PROJECT_NAME} PRIVATE ${SAMRAI_LIBRARIES})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
