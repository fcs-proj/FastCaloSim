# Set up the project.
cmake_minimum_required( VERSION 3.5 )

set(CMAKE_CXX_STANDARD 17)

project( "ATLASMagneticFieldMapPlugin" VERSION 1.0 LANGUAGES CXX )

# Find the header and source files.
file( GLOB SOURCES src/*.cxx )
file(GLOB HEADERS src/*.h)

set(PROJECT_SOURCES ${SOURCES} ${HEADERS})

# Set up the library.
add_library(ATLASMagneticFieldMapPlugin SHARED ${SOURCES})

find_package(Eigen3 REQUIRED)
find_package(Geant4 REQUIRED)

message( STATUS "Found Geant4: ${Geant4_INCLUDE_DIR}")

include(${Geant4_USE_FILE})

# Use the GNU install directory names.
include( GNUInstallDirs )

target_link_libraries ( ATLASMagneticFieldMapPlugin PUBLIC ${CMAKE_DL_LIBS} ${Geant4_LIBRARIES} Eigen3::Eigen)

set_target_properties( ATLASMagneticFieldMapPlugin PROPERTIES
		       VERSION ${PROJECT_VERSION}
		       SOVERSION ${PROJECT_VERSION_MAJOR} )
			   
			
install( TARGETS ATLASMagneticFieldMapPlugin
	 LIBRARY DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MagneticFieldPlugins
	 COMPONENT Runtime
	 NAMELINK_COMPONENT Development )