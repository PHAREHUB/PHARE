
cmake_minimum_required(VERSION 3.0)
project(initializer)




set(SOURCE_INC data_provider.h
               python_data_provider.h
               restart_data_provider.h)


set(SOURCE_CPP python_data_provider.cpp data_provider.cpp)



find_package(pybind11 CONFIG QUIET)

if (NOT pybind11_FOUND)

    message("downloading subproject pybind11")
    set(PYBIND11_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pybind11)

    if (NOT EXISTS ${PYBIND11_SRCDIR})
      execute_process(
        COMMAND ${Git} clone https://github.com/pybind/pybind11 ${PYBIND11_SRCDIR}
        )
    endif()


        add_subdirectory(${PYBIND11_SRCDIR})

endif()

set(PYBIND11_CPP_STANDARD -std=c++17)

add_library(initializer SHARED ${SOURCE_CPP} ${SOURCE_INC})
target_compile_options(initializer PRIVATE ${PHARE_WERROR_FLAGS})



target_include_directories(${PROJECT_NAME}  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../subprojects>
  $<INSTALL_INTERFACE:include/phare/initializer>)

target_link_libraries(initializer PUBLIC pybind11::embed)


# this is the python module that wraps the add() method
pybind11_add_module(pyphare pyphare.cpp)
target_link_libraries(pyphare PUBLIC initializer)


