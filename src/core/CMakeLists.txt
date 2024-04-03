cmake_minimum_required (VERSION 3.20.1)
project(phare_core)

set( SOURCES_INC
     data/electromag/electromag.hpp
     data/field/field.hpp
     data/grid/gridlayoutdefs.hpp
     data/grid/gridlayout.hpp
     data/grid/gridlayout_impl.hpp
     data/grid/gridlayoutimplyee.hpp
     data/grid/gridlayout_utils.hpp
     data/ndarray/ndarray_vector.hpp
     data/particles/particle.hpp
     data/particles/particle_utilities.hpp
     data/particles/particle_array.hpp
     data/ions/ion_population/particle_pack.hpp
     data/ions/ion_population/ion_population.hpp
     data/ions/ions.hpp
     data/electrons/electrons.hpp
     data/ions/particle_initializers/particle_initializer.hpp
     data/ions/particle_initializers/maxwellian_particle_initializer.hpp
     data/ions/particle_initializers/particle_initializer_factory.hpp
     data/tensorfield/tensorfield.hpp
     data/vecfield/vecfield.hpp
     data/vecfield/vecfield_component.hpp
     data/vecfield/vecfield_initializer.hpp
     hybrid/hybrid_quantities.hpp
     pic/pic_quantities.hpp
     numerics/boundary_condition/boundary_condition.hpp
     numerics/interpolator/interpolator.hpp
     numerics/pusher/boris.hpp
     numerics/pusher/pusher.hpp
     numerics/pusher/pusher_factory.hpp
     numerics/ampere/ampere.hpp
     numerics/faraday/faraday.hpp
     numerics/ohm/ohm.hpp
     numerics/moments/moments.hpp
     numerics/ion_updater/ion_updater.hpp
     models/physical_state.hpp
     models/hybrid_state.hpp
     models/mhd_state.hpp
     models/pic_state.hpp
     utilities/box/box.hpp
     utilities/algorithm.hpp
     utilities/constants.hpp
     utilities/index/index.hpp
     utilities/meta/meta_utilities.hpp
     utilities/partitionner/partitionner.hpp
     utilities/point/point.hpp
     utilities/range/range.hpp
     utilities/types.hpp
     utilities/mpi_utils.hpp
   )

set( SOURCES_CPP
     data/ions/particle_initializers/maxwellian_particle_initializer.cpp
     utilities/index/index.cpp
     utilities/mpi_utils.cpp
    )

find_package(MPI)

add_library(${PROJECT_NAME}  ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME}  PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME}  PRIVATE phare_initializer ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    PUBLIC ${PHARE_BASE_LIBS}
  )
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
target_include_directories(${PROJECT_NAME}  PUBLIC ${MPI_C_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>)
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

