cmake_minimum_required (VERSION 3.9)

project(phare_bench_solver_ppc)

add_phare_cpp_benchmark(11 ${PROJECT_NAME} solver_ppc_bench ${CMAKE_CURRENT_BINARY_DIR} solver_ppc_bench.h)
