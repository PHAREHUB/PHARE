cmake_minimum_required (VERSION 3.9)

project(phare_bench_solver_ppc)

add_phare_cpp_benchmark(11 ${PROJECT_NAME} ion_updater_bench ${CMAKE_CURRENT_BINARY_DIR})

#add_phare_cpp_benchmark(11 ${PROJECT_NAME} ion_updater_bench_omp ${CMAKE_CURRENT_BINARY_DIR})
