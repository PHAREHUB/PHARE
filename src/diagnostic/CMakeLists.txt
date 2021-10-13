cmake_minimum_required (VERSION 3.9)

project(phare_diagnostic)

set(SOURCES_INC
  ${PROJECT_SOURCE_DIR}/diagnostic_manager.h
  ${PROJECT_SOURCE_DIR}/diagnostic_writer.h
  ${PROJECT_SOURCE_DIR}/diagnostic_props.h
)

if (HighFive)
 set(SOURCES_INC
   ${SOURCES_INC}
   ${PROJECT_SOURCE_DIR}/detail/h5writer.h
   ${PROJECT_SOURCE_DIR}/detail/h5_utils.h
   ${PROJECT_SOURCE_DIR}/detail/h5typewriter.h
   ${PROJECT_SOURCE_DIR}/detail/types/particle.h
   ${PROJECT_SOURCE_DIR}/detail/types/electromag.h
   ${PROJECT_SOURCE_DIR}/detail/types/fluid.h
   ${PROJECT_SOURCE_DIR}/detail/types/meta.h
 )
endif()


add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME}  INTERFACE
  $<INSTALL_INTERFACE:include/phare/diagnostic>)

target_link_libraries(${PROJECT_NAME} INTERFACE phare_amr)

