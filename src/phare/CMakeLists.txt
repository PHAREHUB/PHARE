
cmake_minimum_required(VERSION 3.0)
project(phare-exe)

add_executable(phare-exe phare.cpp)
target_compile_options(phare-exe PRIVATE ${PHARE_WERROR_FLAGS})
target_include_directories(initializer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
)

target_link_libraries(phare-exe PUBLIC phare_core phare_amr initializer phare_diagnostic)
