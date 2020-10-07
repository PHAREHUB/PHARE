cmake_minimum_required (VERSION 3.9)

project(phare_bench_particles)

add_phare_cpp_benchmark(11 ${PROJECT_NAME} interop ${CMAKE_CURRENT_BINARY_DIR})
