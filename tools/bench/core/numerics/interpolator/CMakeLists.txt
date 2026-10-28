cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_interpolator)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} bench_main ${CMAKE_CURRENT_BINARY_DIR})
add_phare_cpp_benchmark(10 ${PROJECT_NAME} interpolator_access ${CMAKE_CURRENT_BINARY_DIR})
