if(NOT TARGET spatialindex)
  find_package(libspatialindex CONFIG QUIET NO_MODULE)

  if(NOT libspatialindex_FOUND)
    message(STATUS "libspatialindex not found, building from source")

    # Only override CXX flags once
    set(_saved_cxx_flags "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${_saved_cxx_flags} -w -O1 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")

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

    # Only populate once
    FetchContent_GetProperties(libspatialindex)

    if(NOT libspatialindex_POPULATED)
      FetchContent_Populate(libspatialindex)
      add_subdirectory(${libspatialindex_SOURCE_DIR} ${libspatialindex_BINARY_DIR})
    endif()

    set(libspatialindex_INCLUDE_DIRS "${libspatialindex_SOURCE_DIR}/include")
    set(libspatialindex_FOUND TRUE)

    set(CMAKE_CXX_FLAGS "${_saved_cxx_flags}")
    set(BUILD_TESTING ON CACHE BOOL "Enable tests again" FORCE)
  endif()
endif()
