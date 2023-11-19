cmake_minimum_required (VERSION 3.9)

project(phare_configurator)

find_package(MPI REQUIRED)
foreach(LIB ${MPI_LIBRARY})
    get_filename_component(MPI_LIBRARY_PATH ${LIB} DIRECTORY)
endforeach()


set (_phare_mpi_libs )

if(DEFINED MPI_CXX_LIB_NAMES)
  set (_phare_mpi_libs ${MPI_CXX_LIB_NAMES})

elseif(DEFINED MPI_C_LIB_NAMES)
  set (_phare_mpi_libs ${MPI_C_LIB_NAMES})

else(DEFINED MPI_C_LIB_NAMES)
  message(ERROR "MPI not found!")

endif()

add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_C_HEADER_DIR})
target_link_libraries(${PROJECT_NAME} ${_phare_mpi_libs})
