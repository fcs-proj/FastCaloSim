# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

cmake_minimum_required(VERSION 3.14)

# Globbing source and header files
file(GLOB g4app_SOURCES "src/*.cc")

# Declare the library
add_library(g4app SHARED ${g4app_SOURCES})

# Make include directory accessible to libraries that link against it
target_include_directories(g4app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

# Link against Geant4 and FastCaloSim (for G4 transpot tool)
target_link_libraries(g4app PUBLIC ${Geant4_LIBRARIES} FastCaloSim::FastCaloSim)
