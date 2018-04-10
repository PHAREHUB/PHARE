cmake_minimum_required(VERSION 3.3)
project(phare_core)

set( SOURCES_INC
     data/field/field.h
     data/vecfield/vecfield.h
     data/ndarray/ndarray_vector.h
     data/particles/particle_array.h
     data/species/particle_pack.h
     data/particles/particle.h
     data/grid/gridlayout.h
     data/electromag/electromag.h
     numerics/interpolator/interpolator.h
     hybrid/hybrid_quantities.h
   )

 set( SOURCES_CPP
      data/ndarray/ndarray_vector.cpp
     # numerics/interpolator/interpolator.cpp
    )

add_library(phare_core ${SOURCES_INC} ${SOURCES_CPP})
target_include_directories(phare_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/core>)

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)

