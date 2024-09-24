# Findlwtnn.cmake - Locate lwtnn library
# Defines:
#  LWTNN_FOUND
#  LWTNN_INCLUDE_DIRS
#  LWTNN_LIBRARIES
#  Imported target lwtnn::lwtnn

unset(LWTNN_INCLUDE_DIR CACHE)
unset(LWTNN_LIBRARY CACHE)
unset(LWTNN_INCLUDE_DIRS CACHE)
unset(LWTNN_LIBRARIES CACHE)
unset(LWTNN_FOUND CACHE)

set(LWTNN_HINTS)

# Check for environment variables
if(DEFINED ENV{LWTNN_HOME})
    list(APPEND LWTNN_HINTS $ENV{LWTNN_HOME})
endif()
if(DEFINED ENV{LWTNN__HOME})
    list(APPEND LWTNN_HINTS $ENV{LWTNN__HOME})
endif()

# Check for CMake variables
if(DEFINED LWTNN_HOME)
    list(APPEND LWTNN_HINTS ${LWTNN_HOME})
endif()
if(DEFINED LWTNN__HOME)
    list(APPEND LWTNN_HINTS ${LWTNN__HOME})
endif()

# Append standard installation prefixes
list(APPEND LWTNN_HINTS /usr /usr/local /opt)

# Find include directory
find_path(LWTNN_INCLUDE_DIR
  NAMES lwtnn/LightweightNeuralNetwork.hh
  HINTS ${LWTNN_HINTS}
  PATH_SUFFIXES include
  DOC "Path to lwtnn include directory"
)

# Find library
find_library(LWTNN_LIBRARY
  NAMES lwtnn
  HINTS ${LWTNN_HINTS}
  PATH_SUFFIXES lib
  DOC "Path to lwtnn library"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(lwtnn
  REQUIRED_VARS LWTNN_INCLUDE_DIR LWTNN_LIBRARY
  VERSION_VAR LWTNN_VERSION
)

if(LWTNN_FOUND)
  set(LWTNN_INCLUDE_DIRS ${LWTNN_INCLUDE_DIR})
  set(LWTNN_LIBRARIES ${LWTNN_LIBRARY})

  # Determine the library type (optional)
  get_filename_component(LWTNN_LIB_EXTENSION "${LWTNN_LIBRARY}" EXT)
  if(LWTNN_LIB_EXTENSION STREQUAL ".a")
    set(LWTNN_LIB_TYPE STATIC)
  else()
    set(LWTNN_LIB_TYPE SHARED)
  endif()

  # Create the imported target
  add_library(lwtnn::lwtnn ${LWTNN_LIB_TYPE} IMPORTED)

  # Set the properties of the imported target
  set_target_properties(lwtnn::lwtnn PROPERTIES
    IMPORTED_LOCATION "${LWTNN_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LWTNN_INCLUDE_DIRS}"
  )
endif()

mark_as_advanced(LWTNN_INCLUDE_DIR LWTNN_LIBRARY LWTNN_INCLUDE_DIRS LWTNN_LIBRARIES)
