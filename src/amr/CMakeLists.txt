cmake_minimum_required(VERSION 3.3)
project(phare_amr)

set( SOURCES_INC
     data/particles/refine/particles_data_split.h
     data/particles/refine/split.h
     data/particles/particles_data.h
     data/particles/particles_data_factory.h
     data/particles/particles_variable.h
     data/field/coarsening/field_coarsen_operator.h
     data/field/coarsening/field_coarsen_index_weight.h
     data/field/coarsening/coarsen_weighter.h
     data/field/coarsening/field_coarsener.h
     data/field/field_data.h
     data/field/field_data_factory.h
     data/field/field_geometry.h
     data/field/field_overlap.h
     data/field/field_variable.h
     data/field/refine/field_linear_refine.h
     data/field/refine/field_refiner.h
     data/field/refine/linear_weighter.h
     data/field/refine/field_refine_operator.h
     data/field/time_interpolate/field_linear_time_interpolate.h
     tools/field_resource.h
     tools/particle_resource.h
     tools/amr_utils.h
     tools/resources_manager.h
     tools/resources_manager_utilities.h
     tools/resources_guards.h
     evolution/integrator/multiphysics_integrator.h
     evolution/solvers/solver.h
     evolution/solvers/solver_ppc.h
     evolution/solvers/solver_mhd.h
     evolution/transactions/transaction.h
     evolution/transactions/hybrid_transaction.h
     evolution/transactions/hybrid_transaction_strategy.h
     evolution/transactions/hybrid_hybrid_transaction_strategy.h
     evolution/transactions/mhd_hybrid_transaction_strategy.h
     evolution/transactions/mhd_transaction.h
     evolution/transactions/transaction_manager.h
     evolution/transactions/transaction_factory.h
     evolution/transactions/transaction_initializer.h
     evolution/transactions/transaction_info.h
     evolution/transactions/hybrid_transaction_info.h
     evolution/transactions/mhd_transaction_info.h
     physical_models/physical_model.h
     physical_models/hybrid_model.h
     physical_models/mhd_model.h
   )
set( SOURCES_CPP
     data/field/refine/linear_weighter.cpp
     data/particles/refine/split.cpp
     tools/amr_utils.cpp
     evolution/solvers/solver.cpp
     evolution/transactions/transaction.cpp
     data/field/coarsening/field_coarsen.cpp
     tools/resources_manager.cpp
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/amr>)


target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_core
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_geom
  SAMRAI_hier
  SAMRAI_math
  SAMRAI_mesh
  SAMRAI_pdat
  SAMRAI_solv
  SAMRAI_tbox
  SAMRAI_xfer
  )

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
