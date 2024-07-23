set(FMT_HEADERS "")
include(FetchContent)
FetchContent_Declare(
        fmt
        GIT_REPOSITORY  https://github.com/fmtlib/fmt.git
        GIT_TAG         11.0.1
        GIT_PROGRESS    TRUE
)
FetchContent_MakeAvailable(fmt)
