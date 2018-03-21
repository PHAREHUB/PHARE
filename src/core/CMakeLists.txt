cmake_minimum_required(VERSION 3.2)
project(phare_core)

set( SOURCES_INC
     data/field/field.h
   )

 set( SOURCES_CPP
      data/field/field.cpp
    )

add_library(phare_core ${SOURCES_INC} ${SOURCES_CPP})
target_include_directories(phare_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/core>)

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)

