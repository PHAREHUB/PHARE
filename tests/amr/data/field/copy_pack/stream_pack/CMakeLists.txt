cmake_minimum_required (VERSION 3.3)

project(test-field-data-stream-pack)


set(SOURCES_INC
  test_stream_pack_centered_ex.h
  test_stream_pack_centered_ey.h
   )

set(SOURCES_CPP
  test_main.cpp
#  test_stream_pack_centered_ex.cpp disabled cause it's b0rked
#  test_stream_pack_centered_ey.cpp disabled cause it's b0rked
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE
  test-field-data
  )


add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
