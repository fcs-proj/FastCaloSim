# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

cmake_minimum_required(VERSION 3.14)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

include(cmake/prelude.cmake)
include(cmake/silencer.cmake)

project(
    FastCaloSim
    VERSION 0.1.0
    DESCRIPTION "An experiment-independent library for fast calorimeter simulation"
    HOMEPAGE_URL "https://fcs-proj.github.io/FastCaloSim/"
    LANGUAGES CXX
)

# Append the path for CMake to look for any custom modules, e.g. Findlwtnn.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Project dependencies ----

# ROOT is a required dependency for FastCaloSim
find_package(ROOT 6.28 REQUIRED)
# Geant4 (with gdml activated) is a required dependency for FastCaloSim
find_package(Geant4 10.6 REQUIRED gdml)
# lwtnn for neural network inference
find_package(lwtnn 2.13 REQUIRED)
# OnnxRuntime for neural network inference
find_package(ONNXRuntime 1.15 REQUIRED)
# LibXml2 for reading XML files
find_package(libxml REQUIRED)
# Fmt for formatting
find_package(fmt REQUIRED)
# CLHEP for random number generation
find_package(CLHEP REQUIRED)
# Eigen for extrapolation
find_package(Eigen3 REQUIRED)
# TBB for Geant4 transport
find_package(TBB REQUIRED)
# libspatialindex for geometry
find_package(libspatialindex REQUIRED)

# Boost for geometry (headers only)
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost not found")
endif()


# ---- Declare library ----

# Create a shared library
add_library(FastCaloSim_FastCaloSim SHARED "")

# Add the source files
add_subdirectory(source/Core)
add_subdirectory(source/Geometry)
add_subdirectory(source/Transport)
add_subdirectory(source/Extrapolation)

add_library(FastCaloSim::FastCaloSim ALIAS FastCaloSim_FastCaloSim)

# Link against external libraries
target_link_libraries(
    FastCaloSim_FastCaloSim
    PUBLIC
        ROOT::Core ROOT::Tree ROOT::Gpad ROOT::Hist ROOT::Physics ROOT::MathCore ROOT::MathMore ROOT::ROOTDataFrame
        fmt::fmt-header-only
        CLHEP::CLHEP
        ${Geant4_LIBRARIES}
        TBB::tbb
        libspatialindex::spatialindex
    PRIVATE
        ONNXRuntime::ONNXRuntime
        ${LWTNN_LIBRARIES}
        ${LibXml2_LIBRARIES}
)

# Add the include directory to the target
target_include_directories(FastCaloSim_FastCaloSim
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include/FastCaloSim
        ${Boost_INCLUDE_DIRS}
        ${LWTNN_INCLUDE_DIRS}
)
# Generate ROOT dictionary for custom classes when building a shared library
ROOT_GENERATE_DICTIONARY(FastCaloSim_dict MODULE FastCaloSim_FastCaloSim LINKDEF include/FastCaloSim/LinkDef.h)

# Deactivate compiler warnings and static analyzers for third party libraries
# TODO: fix FCS warnings and remove FastCaloSim_FastCaloSim from list
deactivate_checks(FastCaloSim_dict FastCaloSim_FastCaloSim fmt::fmt)

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
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_include_directories(
    FastCaloSim_FastCaloSim SYSTEM
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/export>"
)

target_compile_features(FastCaloSim_FastCaloSim PUBLIC cxx_std_17)


# Define the list of targets to install
set(FastCaloSim_TARGETS FastCaloSim_FastCaloSim)

# ---- Parametrization ----

# Create FastCaloSim::Param target
if (FastCaloSim_PARAM_MODE)
    add_subdirectory(param)
endif()

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
    include(cmake/install-rules.cmake)
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
