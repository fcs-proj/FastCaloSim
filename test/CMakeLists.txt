cmake_minimum_required(VERSION 3.14)

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/silencer.cmake)
include(../cmake/tests.cmake)

# ---- Test Data Location ----
add_definitions(-DTEST_PARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/param/")
add_definitions(-DTEST_GEO_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/geometry/")
add_definitions(-DTEST_MAG_FIELD_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/fieldmap/")

# ---- Python Scripts Location ----
add_definitions(-DPYTHON_SCRIPT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/python/")

# ---- Test Output Location ----
add_definitions(-DTEST_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/output/")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output)

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
deactivate_checks(gtest gtest_main 
                  gmock gmock_main
                  FastCaloSim_FastCaloSim
                  FastCaloSim_dict)

# Build a simple geant4 application
# - Used to test the particle transport
# - Uses simplified geometry 
add_subdirectory(g4app)
deactivate_checks(g4app)

# Build the ATLAS magnetic field map plugin
# - Used to transport particles through the simplified geometry
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/ATLASMagneticFieldMapPlugin ${CMAKE_CURRENT_BINARY_DIR}/ATLASMagneticFieldMapPlugin)
deactivate_checks(ATLASMagneticFieldMapPlugin)
# Pass the path to the magnetic field shared library (ATLASMagneticFieldMapPlugin.so) to the g4app
target_compile_definitions(g4app PRIVATE ATLASMagneticFieldMapPluginLib="$<TARGET_FILE:ATLASMagneticFieldMapPlugin>")

# ---- Tests ----

include(GoogleTest)

#Glob test sources
file(GLOB TEST_SOURCES "source/*_tests.cxx")

# Create test executables
foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
  create_test_executable(${TEST_NAME} ${TEST_SOURCE})
endforeach()

#Tests that are Geant4 dependent
target_link_libraries(transport_tests PRIVATE g4app)


# ---- End-of-file commands ----

add_folders(Test)