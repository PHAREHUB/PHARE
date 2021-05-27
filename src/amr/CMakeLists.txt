cmake_minimum_required(VERSION 3.9)
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
     resources_manager/field_resource.h
     resources_manager/particle_resource.h
     resources_manager/amr_utils.h
     resources_manager/resources_manager.h
     resources_manager/resources_manager_utilities.h
     resources_manager/resources_guards.h
     messengers/quantity_communicator.h
     messengers/communicators.h
     messengers/messenger.h
     messengers/hybrid_messenger.h
     messengers/hybrid_messenger_strategy.h
     messengers/hybrid_hybrid_messenger_strategy.h
     messengers/mhd_hybrid_messenger_strategy.h
     messengers/mhd_messenger.h
     messengers/messenger_factory.h
     messengers/GhostBoxGeometryVariableFillPattern.cpp
     messengers/messenger_info.h
     messengers/hybrid_messenger_info.h
     messengers/mhd_messenger_info.h
     types/amr_types.h
     wrappers/hierarchy.h
     wrappers/integrator.h
     tagging/tagger.h
     tagging/tagger_factory.h
     tagging/hybrid_tagger.h
     tagging/hybrid_tagger_strategy.h
     tagging/default_hybrid_tagger_strategy.h
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
set( SOURCES_CPP
     data/field/refine/linear_weighter.cpp
     resources_manager/amr_utils.cpp
     data/field/coarsening/field_coarsen.cpp
     messengers/messenger_factory.cpp
   )

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION ${PHARE_INTERPROCEDURAL_OPTIMIZATION})


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

