
cmake_minimum_required(VERSION 3.0)
project(phare-exe)



add_executable(phare-exe ${SOURCES_INC} phare.cpp)
target_compile_options(phare-exe PRIVATE ${PHARE_WERROR_FLAGS})
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)


target_link_libraries(phare-exe PUBLIC phare_core phare_amr phare_initializer phare_solver phare_simulator)


# this is the script that PHARE runs to prepare pythonpath in case pyphare is
# ran in a virtualenv
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setpythonpath.py.in ${CMAKE_CURRENT_BINARY_DIR}/setpythonpath.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/init.py ${CMAKE_CURRENT_BINARY_DIR}/init.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py ${CMAKE_CURRENT_BINARY_DIR}/job.py @ONLY)
