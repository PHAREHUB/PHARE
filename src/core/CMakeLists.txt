cmake_minimum_required(VERSION 3.3)
project(phare_core)

set( SOURCES_INC
     data/electromag/electromag.h
     data/field/field.h
     data/grid/gridlayoutdefs.h
     data/grid/gridlayout.h
     data/grid/gridlayout_impl.h
     data/grid/gridlayoutimplyee.h
     data/ndarray/ndarray_vector.h
     data/particles/particle.h
     data/particles/particle_array.h
     data/species/particle_pack.h
     data/vecfield/vecfield.h
     data/vecfield/vecfield_component.h
     hybrid/hybrid_quantities.h
     numerics/boundary_condition/boundary_condition.h
     numerics/interpolator/interpolator.h
     numerics/pusher/boris.h
     numerics/pusher/pusher.h
     numerics/pusher/pusher_factory.h
     numerics/ampere/ampere.h
     numerics/faraday/faraday.h
     utilities/box/box.h
     utilities/constants.h
     utilities/index/index.h
     utilities/meta/meta_utilities.h
     utilities/particle_selector/particle_selector.h
     utilities/partitionner/partitionner.h
     utilities/point/point.h
     utilities/range/range.h
     utilities/types.h
   )

set( SOURCES_CPP
     data/ndarray/ndarray_vector.cpp
     utilities/index/index.cpp
    )

add_library(phare_core ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(phare_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/core>)

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)

