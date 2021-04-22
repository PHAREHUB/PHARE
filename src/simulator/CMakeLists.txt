cmake_minimum_required(VERSION 3.9)
project(phare_simulator)

set( SOURCES_INC
    simulator.h
    simulator.cpp
    phare_types.h
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} )
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS} -DPHARE_HAS_HIGHFIVE=${PHARE_HAS_HIGHFIVE})
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/simulator>)


target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
  Threads::Threads
  )

