cmake_minimum_required (VERSION 3.3)

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
 )
endif()

set(SOURCES_CPP
  ${PROJECT_SOURCE_DIR}/diagnostic_manager.cpp
)

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})
target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})

target_include_directories(${PROJECT_NAME}  PUBLIC
  $<INSTALL_INTERFACE:include/phare/diagnostic>)

target_link_libraries(${PROJECT_NAME} PRIVATE phare_initializer phare_solver)

