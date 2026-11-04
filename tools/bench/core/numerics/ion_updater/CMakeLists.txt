cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_ion_updater)

add_phare_cpp_benchmark(10 bench_updater_impl0 bench_ion_updater ${CMAKE_CURRENT_BINARY_DIR})
if(TARGET bench_updater_impl0_bench_ion_updater) # only made in the exec level range
  target_compile_definitions(bench_updater_impl0_bench_ion_updater PRIVATE PHARE_UPDATER_IMPL=0)
endif()

add_phare_cpp_benchmark(10 bench_updater_impl1 bench_ion_updater ${CMAKE_CURRENT_BINARY_DIR})
if(TARGET bench_updater_impl1_bench_ion_updater) # only made in the exec level range
  target_compile_definitions(bench_updater_impl1_bench_ion_updater PRIVATE PHARE_UPDATER_IMPL=1)
endif()
