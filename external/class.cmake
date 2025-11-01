cmake_minimum_required(VERSION 3.13)
include(FetchContent)

FetchContent_Declare(
    class
    GIT_REPOSITORY https://github.com/lesgourg/class_public.git
    GIT_TAG fb762fda9a9b6efce32d41ca0f0a7e99369e7483 # git tag v3.3.3
    GIT_SHALLOW YES
    GIT_PROGRESS TRUE
    USES_TERMINAL_DOWNLOAD TRUE   # <---- this is needed only for Ninja
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/external/CLASS_CMakeLists.txt
        <SOURCE_DIR>/CMakeLists.txt
)

set(FETCHCONTENT_QUIET OFF)
FetchContent_MakeAvailable(class)
