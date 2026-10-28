cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_ndarray)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} bench_ndarray_access ${CMAKE_CURRENT_BINARY_DIR})
