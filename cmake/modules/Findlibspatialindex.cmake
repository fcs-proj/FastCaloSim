if(NOT TARGET spatialindex)
  find_package(libspatialindex CONFIG QUIET)

  if(NOT libspatialindex_FOUND)
    message(STATUS "libspatialindex not found, building from source")

    # Disable tests in the embedded libspatialindex
    set(SIDX_BASE_CPP "" CACHE STRING "" FORCE)
    set(SIDX_COMMON_CXX_FLAGS "" CACHE STRING "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable libspatialindex tests" FORCE)
    include(FetchContent)
    FetchContent_Declare(
      libspatialindex
      GIT_REPOSITORY https://github.com/libspatialindex/libspatialindex
      GIT_TAG 2.1.0
    )

    # Store original flags
    set(_saved_cxx_flags "${CMAKE_CXX_FLAGS}")

    # Override globally *before build*
    set(CMAKE_CXX_FLAGS "${_saved_cxx_flags} -w -O1 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")

    # Ensure position independent code for shared libraries
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    FetchContent_MakeAvailable(libspatialindex)

    target_include_directories(spatialindex PUBLIC
      "$<BUILD_INTERFACE:${libspatialindex_SOURCE_DIR}/include>"
      "$<INSTALL_INTERFACE:include>"
    )
    add_library(libspatialindex::spatialindex ALIAS spatialindex)

    deactivate_checks(spatialindex spatialindex_c)
    # Restore flags
    set(CMAKE_CXX_FLAGS "${_saved_cxx_flags}")

    set(libspatialindex_INCLUDE_DIRS "${libspatialindex_SOURCE_DIR}/include")
    set(libspatialindex_FOUND TRUE)

    set(BUILD_TESTING ON CACHE BOOL "Enable tests again" FORCE)
  endif()
endif()
