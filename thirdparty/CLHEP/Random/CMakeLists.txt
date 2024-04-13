cmake_minimum_required(VERSION 3.14)

# Name of the library
project(CLHEP LANGUAGES CXX)

# Globbing source and header files
file(GLOB CLHEP_SOURCES "*.cxx")

# Declare the library
add_library(CLHEP SHARED ${CLHEP_SOURCES})

# Set the target propertise
set_target_properties(
    CLHEP PROPERTIES
    EXPORT_NAME CLHEP
    OUTPUT_NAME CLHEP
)

# Set the include properties
target_include_directories(
    CLHEP
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/thirdparty>"
)

# Link against ROOT::Core
target_link_libraries(CLHEP PUBLIC ROOT::Core ROOT::MathCore)

# Add CLHEP to the export set
install( TARGETS CLHEP EXPORT FastCaloSimTargets)