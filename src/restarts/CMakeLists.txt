cmake_minimum_required (VERSION 3.20.1)

project(phare_restarts)

set(SOURCES_INC
  ${PROJECT_SOURCE_DIR}/restarts_manager.hpp
  ${PROJECT_SOURCE_DIR}/restarts_writer.hpp
  ${PROJECT_SOURCE_DIR}/restarts_props.hpp
)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE phare_amr)
