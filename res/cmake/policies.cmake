

if(POLICY CMP0074) # hides warning about ${PACKAGE}_ROOT variables
    cmake_policy(SET CMP0074 NEW)
endif()

if (POLICY CMP0069) # allow LTO if available, requires cmake 3.9 (released August 2017)
  cmake_policy(SET CMP0069 NEW)
endif()

