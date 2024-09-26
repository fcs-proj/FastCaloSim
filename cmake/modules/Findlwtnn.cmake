# Locate the lwtnn external package.
#
# Defines:
#  LWTNN_FOUND
#  LWTNN_INCLUDE_DIR
#  LWTNN_INCLUDE_DIRS
#  LWTNN_LIBRARIES
#  LWTNN_LIBRARY_DIRS
#
# The user can set LWTNN_ROOT to guide the script.
#

# Allow user to specify LWTNN_ROOT
if(DEFINED LWTNN_ROOT)
   set(LWTNN_INCLUDE_DIR "${LWTNN_ROOT}/include")
   set(LWTNN_LIBRARY_DIR "${LWTNN_ROOT}/lib")
else()
   # Search for the include directory and library
   find_path(LWTNN_INCLUDE_DIR NAMES "lwtnn" PATH_SUFFIXES include)
   find_library(LWTNN_LIBRARY NAMES lwtnn)
endif()

# If the include directory or library wasn't found, try again with hints
if(NOT LWTNN_INCLUDE_DIR)
   find_path(LWTNN_INCLUDE_DIR NAMES "lwtnn" PATH_SUFFIXES include HINTS ${LWTNN_ROOT})
endif()

if(NOT LWTNN_LIBRARY)
   find_library(LWTNN_LIBRARY NAMES lwtnn HINTS ${LWTNN_LIBRARY_DIR} ${LWTNN_ROOT}/lib)
endif()

# Set the include directories and libraries
set(LWTNN_INCLUDE_DIRS ${LWTNN_INCLUDE_DIR})
set(LWTNN_LIBRARIES ${LWTNN_LIBRARY})

# Handle the standard find_package arguments:
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(lwtnn DEFAULT_MSG LWTNN_INCLUDE_DIRS LWTNN_LIBRARIES)
mark_as_advanced(LWTNN_FOUND LWTNN_INCLUDE_DIR LWTNN_INCLUDE_DIRS LWTNN_LIBRARIES LWTNN_LIBRARY_DIR)
