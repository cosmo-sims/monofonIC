cmake_minimum_required(VERSION 3.13)
include(FetchContent)

FetchContent_Declare(
    class
    GIT_REPOSITORY https://github.com/lesgourg/class_public.git
    GIT_SHALLOW YES
    GIT_PROGRESS TRUE
    USES_TERMINAL_DOWNLOAD TRUE   # <---- this is needed only for Ninja
)

FetchContent_GetProperties(class)
if(NOT class_POPULATED)
    set(FETCHCONTENT_QUIET OFF)
    FetchContent_Populate(class)
    #configure_file(external/CLASS_CMakeLists.txt ${class_SOURCE_DIR}/CMakeLists.txt COPYONLY)
    file(COPY_FILE external/CLASS_CMakeLists.txt ${class_SOURCE_DIR}/CMakeLists.txt )
    add_subdirectory(${class_SOURCE_DIR} ${class_BINARY_DIR})
endif()
