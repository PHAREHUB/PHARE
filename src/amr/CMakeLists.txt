cmake_minimum_required (VERSION 3.20.1)
project(phare_amr)

set( SOURCES_INC
     data/particles/refine/particles_data_split.hpp
     data/particles/refine/split.hpp
     data/particles/particles_data.hpp
     data/particles/particles_data_factory.hpp
     data/particles/particles_variable.hpp
     data/field/coarsening/field_coarsen_operator.hpp
     data/field/coarsening/field_coarsen_index_weight.hpp
     data/field/coarsening/coarsen_weighter.hpp
     data/field/coarsening/default_field_coarsener.hpp
     data/field/coarsening/magnetic_field_coarsener.hpp
     data/field/field_data.hpp
     data/field/field_data_factory.hpp
     data/field/field_geometry.hpp
     data/field/field_overlap.hpp
     data/field/field_variable.hpp
     data/field/refine/field_linear_refine.hpp
     data/field/refine/field_refiner.hpp
     data/field/refine/magnetic_field_refiner.hpp
     data/field/refine/electric_field_refiner.hpp
     data/field/refine/linear_weighter.hpp
     data/field/refine/field_refine_operator.hpp
     data/field/time_interpolate/field_linear_time_interpolate.hpp
     resources_manager/field_resource.hpp
     resources_manager/particle_resource.hpp
     resources_manager/amr_utils.hpp
     resources_manager/resources_manager.hpp
     resources_manager/resources_manager_utilities.hpp
     resources_manager/resources_guards.hpp
     messengers/communicator.hpp
     messengers/refiner_pool.hpp
     messengers/refiner.hpp
     messengers/synchronizer_pool.hpp
     messengers/synchronizer.hpp
     messengers/messenger.hpp
     messengers/hybrid_messenger.hpp
     messengers/hybrid_messenger_strategy.hpp
     messengers/hybrid_hybrid_messenger_strategy.hpp
     messengers/mhd_hybrid_messenger_strategy.hpp
     messengers/mhd_messenger.hpp
     messengers/messenger_factory.hpp
     messengers/messenger_info.hpp
     messengers/hybrid_messenger_info.hpp
     messengers/mhd_messenger_info.hpp
     types/amr_types.hpp
     wrappers/hierarchy.hpp
     wrappers/integrator.hpp
     tagging/tagger.hpp
     tagging/tagger_factory.hpp
     tagging/hybrid_tagger.hpp
     tagging/hybrid_tagger_strategy.hpp
     tagging/default_hybrid_tagger_strategy.hpp
     solvers/solver.hpp
     solvers/solver_ppc.hpp
     solvers/solver_mhd.hpp
     physical_models/physical_model.hpp
     physical_models/hybrid_model.hpp
     physical_models/mhd_model.hpp
     physical_models/pic_model.hpp
     multiphysics_integrator.hpp
     messenger_registration.hpp
     level_initializer/level_initializer.hpp
     level_initializer/hybrid_level_initializer.hpp
     level_initializer/level_initializer_factory.hpp
     data/field/field_variable_fill_pattern.hpp
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

target_link_directories(${PROJECT_NAME} PUBLIC ${MPI_LIBRARY_PATH})
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
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

