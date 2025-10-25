# CMake script to get CLASS git repository version at compile time
# Adapted from monofonIC's version.cmake
# itself taken from Matt Keeter's blog and slightly adapted
# https://www.mattkeeter.com/blog/2018-01-06-versioning/

execute_process(COMMAND git log --pretty=format:'%h' -n 1
                OUTPUT_VARIABLE CLASS_GIT_REV
                WORKING_DIRECTORY ${CLASS_SOURCE_DIR}
                ERROR_QUIET)

# Check whether we got any revision
if ("${CLASS_GIT_REV}" STREQUAL "")
    set(CLASS_GIT_REV "N/A")
    set(CLASS_GIT_DIFF "")
    set(CLASS_GIT_TAG "N/A")
    set(CLASS_GIT_BRANCH "N/A")
else()
    execute_process(
        COMMAND bash -c "git diff --quiet --exit-code || echo +"
        WORKING_DIRECTORY ${CLASS_SOURCE_DIR}
        OUTPUT_VARIABLE CLASS_GIT_DIFF)
    execute_process(
        COMMAND git describe
        WORKING_DIRECTORY ${CLASS_SOURCE_DIR}
        OUTPUT_VARIABLE CLASS_GIT_TAG ERROR_QUIET)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CLASS_SOURCE_DIR}
        OUTPUT_VARIABLE CLASS_GIT_BRANCH)

    string(STRIP "${CLASS_GIT_REV}" CLASS_GIT_REV)
    string(SUBSTRING "${CLASS_GIT_REV}" 1 7 CLASS_GIT_REV)
    string(STRIP "${CLASS_GIT_DIFF}" CLASS_GIT_DIFF)
    string(STRIP "${CLASS_GIT_TAG}" CLASS_GIT_TAG)
    string(STRIP "${CLASS_GIT_BRANCH}" CLASS_GIT_BRANCH)
endif()

if("${CLASS_GIT_TAG}" STREQUAL "")
    set(CLASS_GIT_TAG "N/A")
endif()

set(CLASS_VERSION "const char* CLASS_GIT_REV=\"${CLASS_GIT_REV}${CLASS_GIT_DIFF}\";
const char* CLASS_GIT_TAG=\"${CLASS_GIT_TAG}\";
const char* CLASS_GIT_BRANCH=\"${CLASS_GIT_BRANCH}\";")

if(EXISTS ${CLASS_BINARY_DIR}/class_version.cc)
    file(READ ${CLASS_BINARY_DIR}/class_version.cc CLASS_VERSION_)
else()
    set(CLASS_VERSION_ "")
endif()

if (NOT "${CLASS_VERSION}" STREQUAL "${CLASS_VERSION_}")
    file(WRITE ${CLASS_BINARY_DIR}/class_version.cc "${CLASS_VERSION}")
endif()
