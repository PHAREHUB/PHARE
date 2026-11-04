cmake_minimum_required (VERSION 3.20.1)

project(phare_bench_pusher)

add_phare_cpp_benchmark(10 ${PROJECT_NAME} pusher ${CMAKE_CURRENT_BINARY_DIR})

add_executable(push_raw_particle_gen push_raw_gen.cpp)
target_include_directories(push_raw_particle_gen PRIVATE ${PHARE_PROJECT_DIR})
target_link_libraries(push_raw_particle_gen PUBLIC ${GBENCH_LIBS} phare_simulator)

add_custom_command(
    TARGET push_raw_particle_gen
    COMMAND push_raw_particle_gen
    DEPENDS pusher
)

add_phare_cpp_benchmark(10 push_raw_use_impl0 push_raw_use ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(push_raw_use_impl0_push_raw_use PRIVATE PHARE_UPDATER_IMPL=0)

add_phare_cpp_benchmark(10 push_raw_use_impl1 push_raw_use ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(push_raw_use_impl1_push_raw_use PRIVATE PHARE_UPDATER_IMPL=1)
