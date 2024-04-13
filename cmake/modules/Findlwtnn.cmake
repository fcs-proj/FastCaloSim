# Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
#
# - Locate lwtnn library
# Defines:
#
#  LWTNN_FOUND
#  LWTNN_INCLUDE_DIR
#  LWTNN_INCLUDE_DIRS
#  LWTNN_<component>_LIBRARY
#  LWTNN_<component>_FOUND
#  LWTNN_LIBRARIES
#

set(LWTNN_LIBRARIES)
set(LWTNN_INCLUDE_DIRS)

# Find the lwtnn executable, and set up the binary path using it:
find_package( lwtnn
  PATHS $ENV{LWTNN__HOME}
)

include_directories(${lwtnn_INCLUDE_DIRS})
link_directories(${lwtnn_LIBRARIES})

# Handle the standard find_package arguments:
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( lwtnn
  DEFAULT_MSG
  lwtnn_INCLUDE_DIRS
  lwtnn_LIBRARIES)

mark_as_advanced(LWTNN_FOUND lwtnn_INCLUDE_DIRS
 lwtnn_LIBRARIES)'