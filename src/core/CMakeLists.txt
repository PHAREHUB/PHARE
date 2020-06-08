cmake_minimum_required(VERSION 3.3)
project(phare_core)

set( SOURCES_INC
     data/electromag/electromag.h
     data/field/field.h
     data/grid/gridlayoutdefs.h
     data/grid/gridlayout.h
     data/grid/gridlayout_impl.h
     data/grid/gridlayoutimplyee.h
     data/grid/gridlayout_utils.h
     data/ndarray/ndarray_vector.h
     data/particles/particle.h
     data/particles/particle_utilities.h
     data/particles/particle_array.h
     data/ions/ion_population/particle_pack.h
     data/ions/ion_population/ion_population.h
     data/ions/ions.h
     data/electrons/electrons.h
     data/ions/particle_initializers/particle_initializer.h
     data/ions/particle_initializers/maxwellian_particle_initializer.h
     data/ions/particle_initializers/particle_initializer_factory.h
     data/vecfield/vecfield.h
     data/vecfield/vecfield_component.h
     data/vecfield/vecfield_initializer.h
     hybrid/hybrid_quantities.h
     numerics/boundary_condition/boundary_condition.h
     numerics/interpolator/interpolator.h
     numerics/pusher/boris.h
     numerics/pusher/pusher.h
     numerics/pusher/pusher_factory.h
     numerics/ampere/ampere.h
     numerics/faraday/faraday.h
     numerics/ohm/ohm.h
     numerics/moments/moments.h
     numerics/ion_updater/ion_updater.h
     models/physical_state.h
     models/hybrid_state.h
     models/mhd_state.h
     utilities/box/box.h
     utilities/algorithm.h
     utilities/constants.h
     utilities/index/index.h
     utilities/meta/meta_utilities.h
     utilities/partitionner/partitionner.h
     utilities/point/point.h
     utilities/range/range.h
     utilities/types.h
     utilities/mpi_utils.h
   )

set( SOURCES_CPP
     data/ndarray/ndarray_vector.cpp
     data/ions/particle_initializers/maxwellian_particle_initializer.cpp
     utilities/index/index.cpp
     utilities/mpi_utils.cpp
    )

find_package(MPI)

add_library(${PROJECT_NAME}  ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME}  PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME}  PRIVATE phare_initializer ${MPI_C_LIBRARIES})

target_include_directories(${PROJECT_NAME}  PUBLIC ${MPI_C_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>)


