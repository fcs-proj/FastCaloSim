# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

cmake_minimum_required(VERSION 3.14)

# Globbing source and header files
file(GLOB ATLASFCalGeoPlugin_SOURCES "src/*.cxx")

# Declare the library
add_library(ATLASFCalGeoPlugin SHARED ${ATLASFCalGeoPlugin_SOURCES})

# Make include directory accessible to libraries that link against it
target_include_directories(ATLASFCalGeoPlugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

# Link against FastCaloSim
target_link_libraries(ATLASFCalGeoPlugin PUBLIC FastCaloSim::FastCaloSim)
