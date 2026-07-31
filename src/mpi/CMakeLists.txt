cmake_minimum_required (VERSION 3.20.1)
project(phare_mpi)

set( SOURCES_INC
     mpi_utils.hpp
   )

set( SOURCES_CPP
     mpi_utils.cpp
    )

set( CORE_EXTRA_SOURCES_CPP )
if (withPhlop)
  set( CORE_EXTRA_SOURCES_CPP ${PHARE_PROJECT_DIR}/subprojects/phlop/src/phlop/timing/common_threaded_timer.cpp
   )
endif(withPhlop)

add_library(${PROJECT_NAME}  ${SOURCES_INC} ${SOURCES_CPP} ${CORE_EXTRA_SOURCES_CPP} )
target_compile_options(${PROJECT_NAME}  PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME}  PRIVATE phare_core MPI::MPI_C
    PUBLIC ${PHARE_BASE_LIBS}
  )

# skip mpicxx in compile command if not configured with configurator
if (NOT phare_configurator)
     target_compile_definitions(${PROJECT_NAME} PUBLIC MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX )
endif()

set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")
