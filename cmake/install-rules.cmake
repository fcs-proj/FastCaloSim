if(PROJECT_IS_TOP_LEVEL)
  set(
      CMAKE_INSTALL_INCLUDEDIR "include/FastCaloSim-${PROJECT_VERSION}"
      CACHE STRING ""
  )
  set_property(CACHE CMAKE_INSTALL_INCLUDEDIR PROPERTY TYPE PATH)
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
set(package FastCaloSim)

install(
    DIRECTORY
    include/
    "${PROJECT_BINARY_DIR}/export/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT FastCaloSim_Development
)

install(
    TARGETS FastCaloSim_FastCaloSim
    EXPORT FastCaloSimTargets
    RUNTIME #
    COMPONENT FastCaloSim_Runtime
    LIBRARY #
    COMPONENT FastCaloSim_Runtime
    NAMELINK_COMPONENT FastCaloSim_Development
    ARCHIVE #
    COMPONENT FastCaloSim_Development
    INCLUDES #
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

# Allow package maintainers to freely override the path for the configs
set(
    FastCaloSim_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${package}"
    CACHE STRING "CMake package config location relative to the install prefix"
)
set_property(CACHE FastCaloSim_INSTALL_CMAKEDIR PROPERTY TYPE PATH)
mark_as_advanced(FastCaloSim_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${FastCaloSim_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT FastCaloSim_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${FastCaloSim_INSTALL_CMAKEDIR}"
    COMPONENT FastCaloSim_Development
)

install(
    EXPORT FastCaloSimTargets
    NAMESPACE FastCaloSim::
    DESTINATION "${FastCaloSim_INSTALL_CMAKEDIR}"
    COMPONENT FastCaloSim_Development
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
