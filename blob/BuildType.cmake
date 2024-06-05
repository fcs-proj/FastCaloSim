
# Author: Marcus D. Hanwell
# Source: https://blog.kitware.com/cmake-and-the-default-build-type/

# Set a default build type if none was specified
set(default_build_type "Release")

# TODO: at the moment, we want to build in Release mode by default,
# even if we build from a Git clone, because that is the default mode
# for our users to get the source code.
# But maybe we will want to change this behavior, later?
# if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
#   set(default_build_type "Debug")
# endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  if( COLOR_DEFS )
    message(STATUS "${Blue}INFO: Setting build type to '${default_build_type}' as none was specified.${ColourReset}")
  else()
    message(STATUS "INFO: Setting build type to '${default_build_type}' as none was specified.")
  endif()
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

