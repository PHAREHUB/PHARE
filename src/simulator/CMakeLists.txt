cmake_minimum_required (VERSION 3.20.1)
project(phare_simulator)

set( SOURCES_INC
    simulator.hpp
    simulator.cpp
    phare_types.hpp
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} )
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
)
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")
