cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)
include(cmake/silencer.cmake)


project(
    FastCaloSim
    VERSION 0.1.0
    DESCRIPTION "An experiment-independent library for fast calorimeter simulation"
    HOMEPAGE_URL "https://jbeirer.github.io/FastCaloSim/"
    LANGUAGES CXX
)

# Set the path for CMake to loook for any custom modules, e.g. Findlwtnn.cmake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)  

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Project dependencies ----

# ROOT is a required dependency for FastCaloSim
find_package(ROOT 6.30 CONFIG REQUIRED)
# Geant4 (with gdml activated) is a required dependency for FastCaloSim 
find_package(Geant4 REQUIRED gdml)
# lwtnn for neural network inference
find_package(lwtnn 2.13 CONFIG REQUIRED)
# OnnxRuntime for neural network inference  
find_package(ONNXRuntime 1.16.3 REQUIRED)
# LibXml2 for reading XML files
find_package(PkgConfig REQUIRED)
pkg_check_modules(LibXml2 REQUIRED libxml-2.0)
include_directories(${LibXml2_INCLUDE_DIRS})
link_directories(${LibXml2_LIBRARY_DIRS})
add_definitions(${LibXml2_CFLAGS_OTHER})
# CLHEP for random number generation
find_package(CLHEP REQUIRED)

# HepPDT for particle data table
add_subdirectory(thirdparty/HepPDT)

# ---- Declare library ----

# Create a shared library
add_library(FastCaloSim_FastCaloSim SHARED "")

# Add the source files
add_subdirectory(source/Core)
add_subdirectory(source/Geometry)
add_subdirectory(source/Transport)

add_library(FastCaloSim::FastCaloSim ALIAS FastCaloSim_FastCaloSim)

# Link against external libraries
target_link_libraries(FastCaloSim_FastCaloSim PUBLIC  ROOT::Core ROOT::Tree ROOT::Gpad ROOT::Hist ROOT::Physics ROOT::MathCore ROOT::MathMore 
                                                      CLHEP::CLHEP
                                                      HepPDT 
                                                      lwtnn::lwtnn
                                                      ONNXRuntime::ONNXRuntime 
                                                      ${LibXml2_LIBRARIES}
                                                      ${Geant4_LIBRARIES})

# Add the include directory to the target
target_include_directories(FastCaloSim_FastCaloSim PRIVATE  ${CMAKE_SOURCE_DIR}/include/FastCaloSim)

# Generate ROOT dictionary for custom classes when building a shared library
ROOT_GENERATE_DICTIONARY(FastCaloSim_dict MODULE FastCaloSim_FastCaloSim LINKDEF include/FastCaloSim/LinkDef.h)

# Deactivate compiler warnings and static analyzers for third party libraries
# TODO: fix FCS warnings and remove FastCaloSim_FastCaloSim from list
deactivate_checks(FastCaloSim_dict HepPDT FastCaloSim_FastCaloSim)

include(GenerateExportHeader)
generate_export_header(
    FastCaloSim_FastCaloSim
    BASE_NAME FastCaloSim
    EXPORT_FILE_NAME export/FastCaloSim/FastCaloSim_export.h
    CUSTOM_CONTENT_FROM_VARIABLE pragma_suppress_c4251
)

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
