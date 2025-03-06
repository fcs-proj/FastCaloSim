# Copyright (c) 2025 CERN for the benefit of the FastCaloSimParam project

cmake_minimum_required(VERSION 3.14)

project(
    FastCaloSimParam
    VERSION 0.1.0
    DESCRIPTION "Parameterization for FastCaloSim"
    HOMEPAGE_URL "https://fcs-proj.github.io/FastCaloSim/"
    LANGUAGES CXX
)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/executable.cmake)

if(PROJECT_IS_TOP_LEVEL)
  find_package(FastCaloSim REQUIRED)
endif()

# ---- Declare library ----

# Create a empty shared library
add_library(FastCaloSimParam_FastCaloSimParam SHARED "")

# Add source files
add_subdirectory(source)

add_library(FastCaloSimParam::FastCaloSimParam ALIAS FastCaloSimParam_FastCaloSimParam)

# Link against external libraries
target_link_libraries(
    FastCaloSimParam_FastCaloSimParam
    PUBLIC
        FastCaloSim::FastCaloSim
)

# Generate ROOT dictionary for custom classes when building a shared library
ROOT_GENERATE_DICTIONARY(FastCaloSimParam_dict MODULE FastCaloSimParam_FastCaloSimParam LINKDEF include/FastCaloSimParam/LinkDef.h)

# Deacticate checks for the generated dictionary
deactivate_checks(FastCaloSimParam_dict FastCaloSimParam_FastCaloSimParam)

# Add include directories for the generated dictionary
target_include_directories(
    FastCaloSimParam_FastCaloSimParam ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)


target_compile_features(FastCaloSimParam_FastCaloSimParam PUBLIC cxx_std_17)




# Add executable
add_exec(createParamSlice)

add_folders(param)
