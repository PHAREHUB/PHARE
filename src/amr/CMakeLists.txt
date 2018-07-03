cmake_minimum_required(VERSION 3.3)
project(phare_amr)

set( SOURCES_INC
     data/particles/refine/particles_data_split_on_coarse_boundary.h
     data/particles/refine/split.h
     data/particles/particles_data.h
     data/particles/particles_data_factory.h
     data/particles/particles_variable.h
     data/field/coarsening/field_data_coarsen.h
     data/field/coarsening/field_coarsen.h
     data/field/field_data.h
     data/field/field_data_factory.h
     data/field/field_geometry.h
     data/field/field_overlap.h
     data/field/field_variable.h
     data/field/refine/field_linear_refine.h
     data/field/refine/field_data_linear_refine.h
     tools/field_resource.h
     tools/particle_resource.h
     tools/amr_utils.h
     tools/resources_manager.h
     tools/resources_manager_utilities.h
     tools/resources_guards.h
   )
set( SOURCES_CPP
     data/field/refine/field_linear_refine.cpp
     data/particles/refine/split.cpp
     tools/amr_utils.cpp
     data/field/coarsening/field_coarsen.cpp
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
