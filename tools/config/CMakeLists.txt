cmake_minimum_required (VERSION 3.20.1) # released April 8, 2021 - https://www.kitware.com/cmake-3-20-1-available-for-download/

project(phare_configurator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS C)
find_package(HDF5 REQUIRED COMPONENTS C)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} MPI::MPI_C HDF5::HDF5)
