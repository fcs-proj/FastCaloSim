# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

if( COLOR_DEFS )
  message(STATUS "-----")
  message(STATUS "${BoldYellow}Building with type: ${CMAKE_BUILD_TYPE}${ColourReset}")
  message(STATUS "${BoldYellow}Using C++ standard: ${CMAKE_CXX_STANDARD}${ColourReset}")
  message(STATUS "-----")
else()
  message(STATUS "-----")
  message(STATUS "Building with type: ${CMAKE_BUILD_TYPE}")
  message(STATUS "Using C++ standard: ${CMAKE_CXX_STANDARD}")
  message(STATUS "-----")
endif()
