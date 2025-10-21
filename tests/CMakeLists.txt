# CMake configuration for monofonIC regression tests

# Find Python3 for the comparison script
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Check that h5py is available
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import h5py"
    RESULT_VARIABLE H5PY_CHECK
    OUTPUT_QUIET
    ERROR_QUIET
)

# Check that numpy is available
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

# Set paths
set(TEST_CONFIG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/configs")
set(TEST_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/references")
set(TEST_SCRIPT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
set(COMPARE_SCRIPT "${TEST_SCRIPT_DIR}/compare_hdf5.py")

# Helper function to add a regression test
# Usage: add_regression_test(<test_name> <config_file> <output_file>)
function(add_regression_test TEST_NAME CONFIG_FILE OUTPUT_FILE)
    set(CONFIG_PATH "${TEST_CONFIG_DIR}/${CONFIG_FILE}")
    set(REFERENCE_PATH "${TEST_REFERENCE_DIR}/${OUTPUT_FILE}")
    set(OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_FILE}")

    # Add test that runs monofonIC and compares output
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} -E env
            ${CMAKE_COMMAND}
            -DMONOFONIC_EXECUTABLE=$<TARGET_FILE:monofonIC>
            -DCONFIG_FILE=${CONFIG_PATH}
            -DOUTPUT_FILE=${OUTPUT_PATH}
            -DREFERENCE_FILE=${REFERENCE_PATH}
            -DCOMPARE_SCRIPT=${COMPARE_SCRIPT}
            -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    # Set test properties
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 300  # 5 minute timeout
        LABELS "regression"
    )

    # Check if reference file exists
    if(NOT EXISTS ${REFERENCE_PATH})
        set_tests_properties(${TEST_NAME} PROPERTIES
            DISABLED TRUE
        )
        message(STATUS "Test ${TEST_NAME}: Reference file not found - test will be disabled")
        message(STATUS "  Generate references with: cd build/tests && bash ${TEST_SCRIPT_DIR}/generate_references.sh")
    else()
        message(STATUS "Test ${TEST_NAME}: Registered")
    endif()
endfunction()

# Register the 5 regression tests
add_regression_test(
    test_1lpt_sc_generic
    test_1lpt_sc_generic.conf
    test_1lpt_sc_generic.hdf5
)

add_regression_test(
    test_2lpt_sc_gadget
    test_2lpt_sc_gadget.conf
    test_2lpt_sc_gadget.hdf5
)

add_regression_test(
    test_3lpt_bcc_swift
    test_3lpt_bcc_swift.conf
    test_3lpt_bcc_swift.hdf5
)

add_regression_test(
    test_2lpt_baryons_generic
    test_2lpt_baryons_generic.conf
    test_2lpt_baryons_generic.hdf5
)

add_regression_test(
    test_2lpt_baryons_vrel_gadget
    test_2lpt_baryons_vrel_gadget.conf
    test_2lpt_baryons_vrel_gadget.hdf5
)

# Add a test to verify the comparison script itself works
add_test(
    NAME compare_script_help
    COMMAND ${Python3_EXECUTABLE} ${COMPARE_SCRIPT} --help
)
set_tests_properties(compare_script_help PROPERTIES
    TIMEOUT 10
    LABELS "sanity"
)

#########################################################################################
# MPI consistency tests
#########################################################################################

# Check if MPI is available (same check as main CMakeLists.txt)
if(MPI_CXX_FOUND)
    set(MPI_TEST_SCRIPT "${TEST_SCRIPT_DIR}/test_mpi_consistency.sh")
    set(MPI_TEST_CONFIG "${TEST_CONFIG_DIR}/test_2lpt_mpi.conf")

    # Add MPI consistency test
    add_test(
        NAME test_mpi_consistency
        COMMAND bash ${MPI_TEST_SCRIPT}
                $<TARGET_FILE:monofonIC>
                ${MPI_TEST_CONFIG}
                ${COMPARE_SCRIPT}
                ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    set_tests_properties(test_mpi_consistency PROPERTIES
        TIMEOUT 600  # 10 minute timeout (runs multiple times)
        LABELS "mpi;regression"
    )

    message(STATUS "Test test_mpi_consistency: Registered (tests with 1, 2, and 4 MPI tasks)")
else()
    message(STATUS "Test test_mpi_consistency: Skipped (MPI not available)")
endif()
