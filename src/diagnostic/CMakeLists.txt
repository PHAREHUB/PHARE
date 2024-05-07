cmake_minimum_required (VERSION 3.20.1)

project(phare_diagnostic)

set(SOURCES_INC
  ${PROJECT_SOURCE_DIR}/diagnostic_manager.hpp
  ${PROJECT_SOURCE_DIR}/diagnostic_writer.hpp
  ${PROJECT_SOURCE_DIR}/diagnostic_props.hpp
)

if (HighFive)
 set(SOURCES_INC
   ${SOURCES_INC}
   ${PROJECT_SOURCE_DIR}/detail/h5writer.hpp
   ${PROJECT_SOURCE_DIR}/detail/h5typewriter.hpp
   ${PROJECT_SOURCE_DIR}/detail/types/particle.hpp
   ${PROJECT_SOURCE_DIR}/detail/types/electromag.hpp
   ${PROJECT_SOURCE_DIR}/detail/types/fluid.hpp
   ${PROJECT_SOURCE_DIR}/detail/types/meta.hpp
 )
endif()

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE phare_amr)
