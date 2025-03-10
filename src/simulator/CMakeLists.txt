cmake_minimum_required (VERSION 3.20.1)
project(phare_simulator)

set( SOURCES_INC
    simulator.hpp
    phare_types.hpp
   )

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
)
