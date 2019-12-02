cmake_minimum_required(VERSION 3.3)
project(phare_solver)

set( SOURCES_INC
     solvers/solver.h
     solvers/solver_ppc.h
     solvers/solver_mhd.h
     physical_models/physical_model.h
     physical_models/hybrid_model.h
     physical_models/mhd_model.h
     multiphysics_integrator.h
     messenger_registration.h
   )
set( SOURCES_CPP
     solvers/solver.cpp
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/solver>)


target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_core
  phare_amr
  )


target_include_directories(${PROJECT_NAME}  PUBLIC
  $<INSTALL_INTERFACE:include/phare/amr>)

