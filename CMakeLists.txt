cmake_minimum_required (VERSION 3.3)

project(PHARE VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

set(PHARE_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(test "Build test with google test" ON)
option(coverage "Generate coverage" ON)
option(documentation "Add doxygen target to generate documentation" ON)
option(cppcheck "Enable cppcheck xml report" ON)

option(asan "build with asan support" OFF)
option(ubsan "build with ubsan support" OFF)
option(msan "build with msan support" OFF)

find_program(Git git)

include(CheckCXXCompilerFlag)

find_package (PythonInterp  3.0 REQUIRED)

set(SUBPROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


#*******************************************************************************
#* Cppcheck option
#*******************************************************************************
if (cppcheck)
  find_program(Cppcheck cppcheck)
  if (NOT Cppcheck-NOTFOUND)

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
      COMMAND ${Cppcheck} --enable=all --std=c++11 --language=c++ --xml
      --project=${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json 2> ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
      )

    add_custom_target(cppcheck-xml DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml)
    find_program(Cppcheck-html cppcheck-htmlreport)
    if (NOT Cppcheck-html-NOTFOUND)

      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml/index.html
        DEPENDS cppcheck-xml
        COMMAND ${Cppcheck-html} --file=${CMAKE_CURRENT_BINARY_DIR}/cppcheck.xml
                                 --report-dir=${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml
        )

      add_custom_target(cppcheck-html DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cppcheckHtml/index.html)
    endif()

  endif()
endif()

#*******************************************************************************
#* SAMRAI build option
#*******************************************************************************
find_package(SAMRAI CONFIG QUIET)
if (NOT SAMRAI_FOUND)

  set(SAMRAI_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/samrai)
  set(SAMRAI_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/samrai)

  if (NOT EXISTS ${SAMRAI_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/monwarez/SAMRAI ${SAMRAI_SRCDIR} --branch add-install --recursive
      )
  endif()

  option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable SAMRAI Test so that we can use the googletest pulled after

  add_subdirectory(${SAMRAI_SRCDIR})



endif()



#*******************************************************************************
#* Test and coverage option
#*******************************************************************************
if (test)

  set(GOOGLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/googletest)

  if (NOT EXISTS ${GOOGLE_TEST_DIR})
     execute_process(
     COMMAND ${Git} clone https://github.com/google/googletest ${GOOGLE_TEST_DIR}
     )
  endif()

  add_subdirectory(subprojects/googletest)

  enable_testing()


  IF(coverage)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -Wall -W -Wshadow -Wunused-variable \
          -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
          -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

      add_custom_target(build-time-make-directory ALL
          COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage)

      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --html  --html-details -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          )

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
          COMMAND gcovr --exclude='.*subprojects.*' --exclude='.*tests.*' --exclude='/usr/include/.*' --object-directory ${CMAKE_BINARY_DIR}  -r ${CMAKE_SOURCE_DIR} --xml -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
          )

      add_custom_target(gcovr
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr ${CMAKE_CURRENT_BINARY_DIR}/coverage/coverage.xml
      )

      if(APPLE)
          set(OPPEN_CMD open)
      elseif(UNIX)
          set(OPPEN_CMD xdg-open)
      endif(APPLE)

      add_custom_target(show_coverage
          COMMAND ${OPPEN_CMD} ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html
          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/coverage/index.html gcovr
          )
  ENDIF(coverage)


  add_subdirectory(tests/amr/data/particles)
  add_subdirectory(tests/amr/data/field/coarsening)
  add_subdirectory(tests/amr/data/field/copy_pack)
  add_subdirectory(tests/amr/data/field/geometry)
  add_subdirectory(tests/amr/data/field/overlap)
  add_subdirectory(tests/amr/data/field/refine)
  add_subdirectory(tests/amr/data/field/variable)
  add_subdirectory(tests/amr/data/field/time_interpolate)
  add_subdirectory(tests/amr/tools/resources_manager)
  add_subdirectory(tests/amr/transactions)
  add_subdirectory(tests/amr/models)
  add_subdirectory(tests/amr/multiphysics_integrator)
  add_subdirectory(tests/core/data/ndarray)
  add_subdirectory(tests/core/data/field)
  add_subdirectory(tests/core/data/gridlayout)
  add_subdirectory(tests/core/data/vecfield)
  add_subdirectory(tests/core/data/particles)
  add_subdirectory(tests/core/data/ions)
  add_subdirectory(tests/core/data/ion_population)
  add_subdirectory(tests/core/data/fluid_particle_initializer)
  add_subdirectory(tests/core/utilities/box)
  add_subdirectory(tests/core/utilities/particle_selector)
  add_subdirectory(tests/core/utilities/partitionner)
  add_subdirectory(tests/core/utilities/range)
  add_subdirectory(tests/core/utilities/index)
  add_subdirectory(tests/core/numerics/boundary_condition)
  add_subdirectory(tests/core/numerics/interpolator)
  add_subdirectory(tests/core/numerics/pusher)
  add_subdirectory(tests/core/numerics/ampere)
  add_subdirectory(tests/core/numerics/faraday)

endif()



#*******************************************************************************
#* Here the Core and AMR lib
#*******************************************************************************

add_subdirectory(src/core)
add_subdirectory(src/amr)

#*******************************************************************************
#* Documentation option
#*******************************************************************************
if (documentation)
    find_package(Doxygen)
    option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

    if(BUILD_DOCUMENTATION)
        if(NOT DOXYGEN_FOUND)
            message(FATAL_ERROR "Doxygen is needed to build the documentation.")
        endif()

        set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
        set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # The project version number.
        set(VERSION_MAJOR   1  )#CACHE STRING "Project major version number.")
        set(VERSION_MINOR   0  )#CACHE STRING "Project minor version number.")
        set(VERSION_PATCH   0  )#CACHE STRING "Project patch version number.")
        set(doxy_main_page "doc/phare.md")
        mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

        configure_file(${doxyfile_in} ${doxyfile} @ONLY)

        add_custom_target(doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}  ${CMAKE_SOURCE_DIR}/doc ${CMAKE_SOURCE_DIR}/src  ${CMAKE_SOURCE_DIR}/tests
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)

        install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
    endif()



  set(DOC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/pharead)

  if (NOT EXISTS ${DOC_DIR})
    execute_process(
      COMMAND ${Git} clone https://hephaistos.lpp.polytechnique.fr/rhodecode/GIT_REPOSITORIES/LPP/phare/pharead
                     ${DOC_DIR}
      )
  endif()

  add_subdirectory(subprojects/pharead)


endif(documentation)


