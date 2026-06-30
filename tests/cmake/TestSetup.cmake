find_package(Python3 COMPONENTS Interpreter REQUIRED)

execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import h5py"
    RESULT_VARIABLE H5PY_CHECK
    OUTPUT_QUIET
    ERROR_QUIET
)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE NUMPY_CHECK
    OUTPUT_QUIET
    ERROR_QUIET
)

if(NOT H5PY_CHECK EQUAL 0 OR NOT NUMPY_CHECK EQUAL 0)
    if(NOT H5PY_CHECK EQUAL 0)
        message(WARNING "Python h5py module not found. Tests will be registered but may fail.")
    endif()
    if(NOT NUMPY_CHECK EQUAL 0)
        message(WARNING "Python numpy module not found. Tests will be registered but may fail.")
    endif()
    message(WARNING "Install missing dependencies with: pip3 install h5py numpy")
    message(WARNING "Or with system package manager: brew install py3-h5py py3-numpy (macOS)")
endif()

set(TEST_CONFIG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/configs")
set(TEST_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/references")
set(TEST_SCRIPT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
set(COMPARE_SCRIPT "${TEST_SCRIPT_DIR}/compare_hdf5.py")
