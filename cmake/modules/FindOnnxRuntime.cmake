# Find the ONNX Runtime include directory and library.
#
# This module defines the `onnxruntime` imported target that encodes all
# necessary information in its target properties.


find_library(OnnxRuntime_LIBRARY onnxruntime HINTS /usr/local/lib64/)
  
find_path(OnnxRuntime_INCLUDE_DIR onnxruntime_cxx_api.h HINTS /usr/local/include/onnxruntime)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  OnnxRuntime
  REQUIRED_VARS OnnxRuntime_LIBRARY OnnxRuntime_INCLUDE_DIR)

add_library(OnnxRuntime SHARED IMPORTED)
set_property(TARGET OnnxRuntime PROPERTY IMPORTED_LOCATION ${OnnxRuntime_LIBRARY})
set_property(TARGET OnnxRuntime PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${OnnxRuntime_INCLUDE_DIR})

mark_as_advanced(OnnxRuntime_FOUND OnnxRuntime_INCLUDE_DIR OnnxRuntime_LIBRARY)