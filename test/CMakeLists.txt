cmake_minimum_required(VERSION 3.14)

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/silencer.cmake)
include(../cmake/tests.cmake)

# ---- Test Data ----
add_definitions(-DTEST_PARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/param/")
add_definitions(-DTEST_GEO_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/geometry/")

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

# ---- Tests ----

# Note: 
# - we generate a single test executable which can take advantage of gtest fixture, etc.
# - you can also run with ctest to get a more detailed output, but no fixture support
# - Run a specific test within a test suite as ctest -R TestSuiteName.TestName --verbose

include(GoogleTest)

# Glob test sources
file(GLOB TEST_SOURCES "source/*_tests.cxx")

# Create a test executablecd 
create_test_executable(run_tests "${TEST_SOURCES}")

# ---- End-of-file commands ----

add_folders(Test)