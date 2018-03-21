if (asan)

  set(ASAN_FLAGS "-fsanitize=address" )
  set(CMAKE_REQUIRED_FLAGS ${ASAN_FLAGS})
  check_cxx_compiler_flag( ${ASAN_FLAGS} ADDRESS_SANITIZER)

  if (${ADDRESS_SANITIZER})
    target_compile_options(${PROJECT_NAME} PUBLIC ${ASAN_FLAGS} -fno-omit-frame-pointer)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS ${ASAN_FLAGS})
  else()
    message(FATAL_ERROR "Your compiler: ${CMAKE_CXX_COMPILER_ID} seems to not support asan flags")
  endif()

  unset(CMAKE_REQUIRED_FLAGS)

elseif (ubsan)

  set(UBSAN_FLAGS "-fsanitize=undefined")
  set(CMAKE_REQUIRED_FLAGS ${UBSAN_FLAGS})
  check_cxx_compiler_flag(${UBSAN_FLAGS} UBSAN_SANITIZER)

  if (${UBSAN_SANITIZER})
    target_compile_options(${PROJECT_NAME} PUBLIC ${UBSAN_FLAGS} )
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS ${UBSAN_FLAGS})
  else ()
    message(FATAL_ERROR "Your compiler: ${CMAKE_CXX_COMPILER_ID} seems to not support ubsan flags")
  endif()

  unset(CMAKE_REQUIRED_FLAGS)

elseif(msan)

  set(MSAN_FLAGS "-fsanitize=memory ")
  set(CMAKE_REQUIRED_FLAGS ${MSAN_FLAGS})
  check_cxx_compiler_flag(-fsanitize=memory MEMORY_SANITIZER)

  if (${MEMORY_SANITIZER})
    target_compile_options(${PROJECT_NAME} PUBLIC ${MSAN_FLAGS} -fno-omit-frame-pointer)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS ${MSAN_FLAGS})
  else ()
    message(FATAL_ERROR "Your compiler: ${CMAKE_CXX_COMPILER_ID} seems to not support msan flags")
  endif()

  unset(CMAKE_REQUIRED_FLAGS)

endif()
