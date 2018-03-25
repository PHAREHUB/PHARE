cmake_minimum_required(VERSION 3.2)
project(phare_core)

set( SOURCES_INC
     data/field/field.h
     data/vecfield/vecfield.h
     data/ndarray/ndarray_vector.h
     hybrid/hybrid_quantities.h
   )

 set( SOURCES_CPP
      data/ndarray/ndarray_vector.cpp
    )

add_library(phare_core ${SOURCES_INC} ${SOURCES_CPP})
target_include_directories(phare_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/core>)

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)

