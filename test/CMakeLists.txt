cmake_minimum_required(VERSION 3.14)

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/silencer.cmake)
include(../cmake/tests.cmake)

# ---- Test Data ----
add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")

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

include(GoogleTest)

# Create the test suites. Note test suites need to 
# --> be in the source directory
# --> have the name *_tests.cxx
create_test_suites(param_file truth_state sim_state)

# ---- End-of-file commands ----

add_folders(Test)