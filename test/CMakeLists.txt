cmake_minimum_required(VERSION 3.14)

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/silencer.cmake)
include(../cmake/tests.cmake)

# ---- Test Base Directory ----
add_definitions(-DTEST_BASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

# ---- Create test output Location ----
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output)
add_definitions(-DTEST_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/output/")
# ---- Location to store test references ----
add_definitions(-DTEST_REFS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/refs/")

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

# Boost is used for transport / extrapolation data serialization in the transport test
find_package(Boost REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})

# Link tests against plugins and g4app
target_link_libraries(BasicExtrapolTests PRIVATE ATLASFCalGeoPlugin)
target_link_libraries(BasicSimTests PRIVATE ATLASFCalGeoPlugin)
target_link_libraries(AtlasTransportTests PRIVATE ATLASFCalGeoPlugin g4app ${Boost_LIBRARIES})
target_link_libraries(AtlasSimTests PRIVATE ATLASFCalGeoPlugin g4app)
target_link_libraries(AtlasGeoTests PRIVATE ATLASFCalGeoPlugin)


# ---- End-of-file commands ----

add_folders(Test)
