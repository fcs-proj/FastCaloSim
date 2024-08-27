cmake_minimum_required(VERSION 3.14)

# Name of the library
project(HepPDT LANGUAGES CXX)

# Globbing source and header files
file(GLOB HepPDT_SOURCES "*.cxx")

# Declare the library
add_library(HepPDT SHARED ${HepPDT_SOURCES})

# Set the target properties
set_target_properties(
    HepPDT PROPERTIES
    EXPORT_NAME HepPDT
    OUTPUT_NAME HepPDT
)

# Set the include properties
target_include_directories(
    HepPDT
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/thirdparty>"
)

# Link against ROOT::Core
target_link_libraries(HepPDT PUBLIC ROOT::Core)

# Add HepPDT to the export set
install( TARGETS HepPDT EXPORT FastCaloSimTargets)