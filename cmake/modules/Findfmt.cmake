set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(FMT_HEADERS "")

# Attempt to find fmt locally
find_package(fmt CONFIG QUIET NO_MODULE)

if (NOT fmt_FOUND)
    message(STATUS "Fmt not found, building from source")
    include(FetchContent)
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 11.0.2
        GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(fmt)
endif()
