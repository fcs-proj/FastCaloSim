# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

cmake_minimum_required(VERSION 3.14)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/silencer.cmake)
include(../cmake/tests.cmake)

# ---- Test Base Directory ----
add_compile_definitions(TEST_BASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

# ---- Define param binary directory ----
add_compile_definitions(PARAM_BINARY_DIR="${CMAKE_BINARY_DIR}/param/")

# ---- Create test output Location ----
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output)
add_compile_definitions(TEST_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/output/")

# ---- Location to store test references ----
add_compile_definitions(TEST_REFS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/refs/")

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(FastCaloSim REQUIRED)
  enable_testing()
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG v1.14.0
  GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(googletest)

# Boost for FastCaloSim and serialization
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# libspatialindex for geometry
find_package(libspatialindex REQUIRED)
include_directories(${libspatialindex_INCLUDE_DIRS})

# Deactivate compiler warnings and static analyzers
# Note: test should be linked with gtest_main unless
# you want to write your own main function
deactivate_checks(gtest gtest_main
                  gmock gmock_main
                  FastCaloSim_FastCaloSim
                  FastCaloSim_dict)

# Build a simple geant4 application
# - Used to test the particle transport
# - Uses simplified geometry
add_subdirectory(g4app)
deactivate_checks(g4app)

# Get the nlohmann json library
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz DOWNLOAD_EXTRACT_TIMESTAMP true)
set(NLOHMANN_JSON_SYSTEM_INCLUDE "" CACHE STRING "" FORCE)
FetchContent_MakeAvailable(json)

# Build the ATLAS geometry plugin
add_subdirectory(TestPlugins/ATLAS/ATLASFCalGeoPlugin)
deactivate_checks(ATLASFCalGeoPlugin)

# Build the ATLAS magnetic field map plugin
# - Used to transport particles through the simplified geometry
add_subdirectory(TestPlugins/ATLAS/ATLASMagneticFieldMapPlugin)
deactivate_checks(ATLASMagneticFieldMapPlugin)
# Pass the path to the magnetic field shared library (MagneticFieldMapPlugin.so) to the g4app
target_compile_definitions(g4app PRIVATE ATLASMagneticFieldMapPluginLib="$<TARGET_FILE:ATLASMagneticFieldMapPlugin>")

# ---- Tests ----

include(GoogleTest)

#Glob test sources
file(GLOB TEST_SOURCES "source/*Tests.cxx")

# Create test executables
foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
  create_test_executable(${TEST_NAME} ${TEST_SOURCE})
endforeach()

# Link tests against plugins and g4app
target_link_libraries(BasicExtrapolTests PRIVATE ATLASFCalGeoPlugin)
target_link_libraries(BasicSimTests PRIVATE ATLASFCalGeoPlugin)
target_link_libraries(AtlasTransportTests PRIVATE ATLASFCalGeoPlugin g4app)
target_link_libraries(AtlasSimTests PRIVATE ATLASFCalGeoPlugin g4app)
target_link_libraries(AtlasGeoTests PRIVATE ATLASFCalGeoPlugin)

# Add python tests
find_package(Python3 COMPONENTS Interpreter REQUIRED)
file(GLOB PYTHON_TESTS "python/*tests.py")
foreach(PYTHON_TEST ${PYTHON_TESTS})
  get_filename_component(TEST_NAME ${PYTHON_TEST} NAME_WE)
  add_test(NAME ${TEST_NAME}
           COMMAND ${Python3_EXECUTABLE} ${PYTHON_TEST})
endforeach()

# ---- End-of-file commands ----

add_folders(Test)
