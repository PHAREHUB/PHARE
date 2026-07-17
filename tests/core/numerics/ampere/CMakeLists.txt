cmake_minimum_required (VERSION 3.20.1)

project(test-ampere)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  ${GTEST_LIBS})

add_custom_command(TARGET test-ampere
                   PRE_BUILD
                   COMMAND  "PYTHONPATH=${PHARE_PYTHONPATH}" ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_ampere.py ${CMAKE_CURRENT_BINARY_DIR})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



if(withKokkos)

    add_executable(amper_gpu_kokkos test_gpu_main.cpp)
    add_no_mpi_phare_test(amper_gpu_kokkos ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(amper_gpu_kokkos PRIVATE
      phare_core
      phare_simulator
      ${GTEST_LIBS})
    target_include_directories(amper_gpu_kokkos PRIVATE
      ${GTEST_INCLUDE_DIRS}
      )

endif(withKokkos)

if(withRaja AND withUmpire)

    add_executable(amper_llnl test_gpu_main.cpp)
    add_no_mpi_phare_test(amper_llnl ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(amper_llnl PRIVATE
      phare_core
      phare_simulator
      ${PHARE_BASE_LIBS}
      ${GTEST_LIBS}
    )
    target_include_directories(amper_llnl PRIVATE
      ${GTEST_INCLUDE_DIRS}
      )

endif(withRaja AND withUmpire)
