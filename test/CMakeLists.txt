cmake_minimum_required(VERSION 3.14)

project(FastCaloSimTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(FastCaloSim REQUIRED)
  enable_testing()
endif()

# ---- Tests ----

add_executable(FastCaloSim_test source/FastCaloSim_test.cxx)
target_link_libraries(FastCaloSim_test PRIVATE FastCaloSim::FastCaloSim)
target_compile_features(FastCaloSim_test PRIVATE cxx_std_17)

add_test(NAME FastCaloSim_test COMMAND FastCaloSim_test)

# ---- End-of-file commands ----

add_folders(Test)