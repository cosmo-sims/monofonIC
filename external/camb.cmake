cmake_minimum_required(VERSION 3.13)
include(FetchContent)

# Find Python3 with interpreter and development components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

# Fetch pybind11 for Python-C++ bindings
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.11.1
    GIT_SHALLOW YES
    GIT_PROGRESS TRUE
    USES_TERMINAL_DOWNLOAD TRUE
)

set(FETCHCONTENT_QUIET OFF)
FetchContent_MakeAvailable(pybind11)

# Check if CAMB Python package is installed
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import camb; print(camb.__version__)"
    RESULT_VARIABLE CAMB_IMPORT_RESULT
    OUTPUT_VARIABLE CAMB_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT CAMB_IMPORT_RESULT EQUAL 0)
    message(STATUS "CAMB Python package not found. Attempting to install...")

    # Try with --user first
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install camb --user
        RESULT_VARIABLE CAMB_INSTALL_RESULT
        OUTPUT_VARIABLE CAMB_INSTALL_OUTPUT
        ERROR_VARIABLE CAMB_INSTALL_ERROR
    )

    # If that fails, try with --break-system-packages (for Homebrew Python)
    if(NOT CAMB_INSTALL_RESULT EQUAL 0)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m pip install camb --user --break-system-packages
            RESULT_VARIABLE CAMB_INSTALL_RESULT
            OUTPUT_VARIABLE CAMB_INSTALL_OUTPUT
            ERROR_VARIABLE CAMB_INSTALL_ERROR
        )
    endif()

    if(NOT CAMB_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR
            "Failed to install CAMB Python package. Please install manually:\n"
            "  ${Python3_EXECUTABLE} -m pip install camb --user\n"
            "Or for Homebrew Python:\n"
            "  ${Python3_EXECUTABLE} -m pip install camb --user --break-system-packages\n"
            "Or use a virtual environment:\n"
            "  python3 -m venv venv && source venv/bin/activate && pip install camb\n"
            "Error: ${CAMB_INSTALL_ERROR}")
    endif()

    # Verify installation
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import camb; print(camb.__version__)"
        RESULT_VARIABLE CAMB_VERIFY_RESULT
        OUTPUT_VARIABLE CAMB_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(NOT CAMB_VERIFY_RESULT EQUAL 0)
        message(FATAL_ERROR "CAMB installation verification failed")
    endif()

    message(STATUS "Successfully installed CAMB ${CAMB_VERSION}")
else()
    message(STATUS "Found CAMB Python package version: ${CAMB_VERSION}")
endif()

# Create interface library for CAMB that links Python and pybind11
add_library(camb_interface INTERFACE)
target_link_libraries(camb_interface INTERFACE
    Python3::Python
    pybind11::embed
)
target_include_directories(camb_interface INTERFACE
    ${Python3_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
)
