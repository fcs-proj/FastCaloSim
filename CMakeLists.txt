cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    FastCaloSim
    VERSION 0.1.0
    DESCRIPTION "An experiment-independent library for fast calorimeter simulation"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Project dependencies ----
find_package(ROOT 6.30 CONFIG REQUIRED)

# ---- Declare library ----

# Globbing all .cxx files in the FastCaloSim directory
file(GLOB FastCaloSim_SOURCES "source/*.cxx")

add_library(
    FastCaloSim_FastCaloSim
    ${FastCaloSim_SOURCES}
)
add_library(FastCaloSim::FastCaloSim ALIAS FastCaloSim_FastCaloSim)

# Link against ROOT libraries
target_link_libraries(FastCaloSim_FastCaloSim PUBLIC ROOT::Core)
# Generate ROOT dictionary for custom classes when building a shared library
if (BUILD_SHARED_LIBS)
  target_include_directories(FastCaloSim_FastCaloSim PRIVATE  ${CMAKE_SOURCE_DIR}/include/FastCaloSim)
  ROOT_GENERATE_DICTIONARY(FastCaloSim_dict FastCaloSim.h MODULE FastCaloSim_FastCaloSim LINKDEF include/FastCaloSim/LinkDef.h)
endif()
# Disable compiler and clang-tidy warnings from ROOT
set_target_properties(FastCaloSim_dict PROPERTIES COMPILE_FLAGS "-w" CXX_CLANG_TIDY "")

include(GenerateExportHeader)
generate_export_header(
    FastCaloSim_FastCaloSim
    BASE_NAME FastCaloSim
    EXPORT_FILE_NAME export/FastCaloSim/FastCaloSim_export.h
    CUSTOM_CONTENT_FROM_VARIABLE pragma_suppress_c4251
)

if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(FastCaloSim_FastCaloSim PUBLIC FASTCALOSIM_STATIC_DEFINE)
endif()

set_target_properties(
    FastCaloSim_FastCaloSim PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME FastCaloSim
    OUTPUT_NAME FastCaloSim
)

target_include_directories(
    FastCaloSim_FastCaloSim ${warning_guard}
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_include_directories(
    FastCaloSim_FastCaloSim SYSTEM
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/export>"
)

target_compile_features(FastCaloSim_FastCaloSim PUBLIC cxx_std_17)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree." "${FastCaloSim_DEVELOPER_MODE}")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

# ---- Developer mode ----

if(NOT FastCaloSim_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of FastCaloSim"
  )
endif()

include(cmake/dev-mode.cmake)
