cmake_minimum_required(VERSION 3.13)
include(FetchContent)

FetchContent_Declare(
    class
    GIT_REPOSITORY https://github.com/lesgourg/class_public.git
    GIT_SHALLOW YES
    GIT_PROGRESS TRUE
    USES_TERMINAL_DOWNLOAD TRUE   # <---- this is needed only for Ninja
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/external/CLASS_CMakeLists.txt
        <SOURCE_DIR>/CMakeLists.txt
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/external/class_version.cmake
        <SOURCE_DIR>/class_version.cmake
)

set(FETCHCONTENT_QUIET OFF)
FetchContent_MakeAvailable(class)
