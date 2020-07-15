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
     level_initializer/level_initializer.h
     level_initializer/hybrid_level_initializer.h
     level_initializer/level_initializer_factory.h
   )

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/amr>
  $<INSTALL_INTERFACE:include/phare/solver>)


target_link_libraries(${PROJECT_NAME} INTERFACE phare_amr)


