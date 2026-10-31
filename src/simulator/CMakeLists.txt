cmake_minimum_required (VERSION 3.20.1)
project(phare_simulator)

set( SOURCES_INC
    simulator.hpp
    simulator_decl.hpp
    simulator_models.hpp
    simulator_solvers.hpp
    simulator_solver_steps.hpp
    simulator_restarts.hpp
    simulator_integrator.hpp
    simulator_messengers.hpp
    simulator_messenger_borders.hpp
    simulator_diagnostics.hpp
    phare_types.hpp
   )

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
)
