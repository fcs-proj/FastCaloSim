# Fetch the polymorphic_variant repository
FetchContent_Declare(polymorphic_variant
    GIT_REPOSITORY https://github.com/Krzmbrzl/polymorphic_variant
    GIT_TAG v1.3.0
)
FetchContent_MakeAvailable(polymorphic_variant)

# Create an INTERFACE library to set up correct include directories
add_library(polymorphic_variant_interface INTERFACE)

# Set the include directories for the new interface target
target_include_directories(polymorphic_variant_interface INTERFACE
    "$<BUILD_INTERFACE:${polymorphic_variant_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

# Install the new interface target
install(TARGETS polymorphic_variant_interface
    EXPORT FastCaloSimTargets
    INCLUDES DESTINATION include
)

# Install the headers
install(DIRECTORY ${polymorphic_variant_SOURCE_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

add_library(polymorphic_variant ALIAS polymorphic_variant_interface)
