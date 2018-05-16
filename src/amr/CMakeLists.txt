cmake_minimum_required(VERSION 3.3)
project(phare_amr)

set( SOURCES_INC
     data/particles/particles_data.h
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


target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_core
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
  )

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
