
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
