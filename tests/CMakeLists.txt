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

#########################################################################################
# Transfer function regression tests
#########################################################################################

# Build test_transfer_functions executable
# Include all source files except main.cc
file(GLOB TEST_TRANSFER_SOURCES
    ${PROJECT_SOURCE_DIR}/src/*.cc
)
list(FILTER TEST_TRANSFER_SOURCES EXCLUDE REGEX ".*main\\.cc$")

# Get plugin files - only transfer function and random plugins (exclude output plugins)
file(GLOB TEST_TRANSFER_PLUGINS
    ${PROJECT_SOURCE_DIR}/src/plugins/transfer_*.cc
    ${PROJECT_SOURCE_DIR}/src/plugins/random_*.cc
)

# Add the test executable (no version.cc needed since we exclude output plugins)
add_executable(test_transfer_functions
    ${CMAKE_CURRENT_SOURCE_DIR}/test_transfer_functions.cc
    ${TEST_TRANSFER_SOURCES}
    ${TEST_TRANSFER_PLUGINS}
)

# Add PANPHASIA sources if enabled (same as main executable)
if(ENABLE_PANPHASIA)
    if(ENABLE_MPI)
        target_sources(test_transfer_functions PRIVATE
            ${PROJECT_SOURCE_DIR}/external/panphasia/panphasia_routines.f
            ${PROJECT_SOURCE_DIR}/external/panphasia/generic_lecuyer.f90
            ${PROJECT_SOURCE_DIR}/external/panphasia_ho/high_order_panphasia_routines.c
            ${PROJECT_SOURCE_DIR}/external/panphasia_ho/pan_mpi_routines.c
            ${PROJECT_SOURCE_DIR}/external/panphasia_ho/uniform_rand_threefry4x64.c
        )
    endif()
endif()

# C++ standard and compile options are inherited from parent CMakeLists.txt

# Include directories
target_include_directories(test_transfer_functions PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external
)

# Link libraries (same as main executable)
target_link_libraries(test_transfer_functions PRIVATE GSL::gsl)

# MPI and FFTW3 linking (same logic as main executable)
if(MPI_CXX_FOUND)
    if(CODE_PRECISION STREQUAL "FLOAT")
        if(FFTW3_SINGLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    elseif(CODE_PRECISION STREQUAL "DOUBLE")
        if(FFTW3_DOUBLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
        if(FFTW3_LONGDOUBLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    endif()
    target_link_libraries(test_transfer_functions PRIVATE MPI::MPI_CXX)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_MPI")
endif()

# FFTW3 serial/threads
if(CODE_PRECISION STREQUAL "FLOAT")
    if(FFTW3_SINGLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_SERIAL)
elseif(CODE_PRECISION STREQUAL "DOUBLE")
    if(FFTW3_DOUBLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_SERIAL)
elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
    if(FFTW3_LONGDOUBLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_SERIAL)
endif()

# HDF5
if(HDF5_FOUND)
    target_link_libraries(test_transfer_functions PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(test_transfer_functions PRIVATE ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(test_transfer_functions PRIVATE "USE_HDF5")
endif()

# PANPHASIA
if(ENABLE_PANPHASIA)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_PANPHASIA")
endif()

# CLASS
if(ENABLE_CLASS)
    target_link_libraries(test_transfer_functions PRIVATE class::libclass)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_CLASS")
endif()

# Helper function to add a transfer function regression test
function(add_transfer_test PLUGIN_NAME CONFIG_NAME)
    set(TEST_NAME "test_transfer_${PLUGIN_NAME}_${CONFIG_NAME}")
    set(CONFIG_FILE "${TEST_CONFIG_DIR}/transfer/${CONFIG_NAME}.conf")
    set(REFERENCE_FILE "${TEST_REFERENCE_DIR}/transfer/${PLUGIN_NAME}/${CONFIG_NAME}.txt")
    set(COMPARE_TF_SCRIPT "${TEST_SCRIPT_DIR}/compare_transfer_functions.py")

    # Add test
    add_test(
        NAME ${TEST_NAME}
        COMMAND $<TARGET_FILE:test_transfer_functions>
                --test
                ${PLUGIN_NAME}
                ${CONFIG_FILE}
                ${REFERENCE_FILE}
                1e-6  # rtol
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    # Set test properties
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 120  # 2 minute timeout
        LABELS "transfer;regression"
    )

    # Check if reference file exists
    if(NOT EXISTS ${REFERENCE_FILE})
        set_tests_properties(${TEST_NAME} PROPERTIES DISABLED TRUE)
        message(STATUS "Test ${TEST_NAME}: Reference file not found - test will be disabled")
    else()
        message(STATUS "Test ${TEST_NAME}: Registered")
    endif()
endfunction()

# Register transfer function tests
# Test configurations
set(TF_CONFIGS fiducial dark_energy massive_nu low_omega_m high_z)

# CLASS plugin tests (all configs)
if(ENABLE_CLASS)
    foreach(CONFIG ${TF_CONFIGS})
        add_transfer_test("CLASS" ${CONFIG})
    endforeach()
endif()

# Eisenstein plugin tests (all configs)
foreach(CONFIG ${TF_CONFIGS})
    add_transfer_test("eisenstein" ${CONFIG})
endforeach()

message(STATUS "")
message(STATUS "Transfer function regression tests configured")
message(STATUS "  Generate references with: cd build && bash ../tests/scripts/generate_transfer_references.sh")
message(STATUS "  Run tests with: ctest -R test_transfer")
