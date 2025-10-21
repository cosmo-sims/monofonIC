# CMake script to run a single monofonIC regression test
# Expected variables:
#   MONOFONIC_EXECUTABLE - path to monofonIC executable
#   CONFIG_FILE - path to test config file
#   OUTPUT_FILE - path where output will be written
#   REFERENCE_FILE - path to reference HDF5 file
#   COMPARE_SCRIPT - path to compare_hdf5.py
#   PYTHON_EXECUTABLE - path to python interpreter

# Check that all required variables are set
if(NOT DEFINED MONOFONIC_EXECUTABLE)
    message(FATAL_ERROR "MONOFONIC_EXECUTABLE not set")
endif()
if(NOT DEFINED CONFIG_FILE)
    message(FATAL_ERROR "CONFIG_FILE not set")
endif()
if(NOT DEFINED OUTPUT_FILE)
    message(FATAL_ERROR "OUTPUT_FILE not set")
endif()
if(NOT DEFINED REFERENCE_FILE)
    message(FATAL_ERROR "REFERENCE_FILE not set")
endif()
if(NOT DEFINED COMPARE_SCRIPT)
    message(FATAL_ERROR "COMPARE_SCRIPT not set")
endif()
if(NOT DEFINED PYTHON_EXECUTABLE)
    message(FATAL_ERROR "PYTHON_EXECUTABLE not set")
endif()

# Check that reference file exists
if(NOT EXISTS ${REFERENCE_FILE})
    message(FATAL_ERROR "Reference file not found: ${REFERENCE_FILE}")
endif()

# Remove any existing output file
file(REMOVE ${OUTPUT_FILE})

# Run monofonIC
message("Running monofonIC with config: ${CONFIG_FILE}")
execute_process(
    COMMAND ${MONOFONIC_EXECUTABLE} ${CONFIG_FILE}
    RESULT_VARIABLE MONOFONIC_RESULT
    OUTPUT_VARIABLE MONOFONIC_OUTPUT
    ERROR_VARIABLE MONOFONIC_ERROR
    TIMEOUT 180  # 3 minute timeout for monofonIC execution
)

# Check that monofonIC ran successfully
if(NOT MONOFONIC_RESULT EQUAL 0)
    message(FATAL_ERROR "monofonIC failed with exit code ${MONOFONIC_RESULT}\nOutput: ${MONOFONIC_OUTPUT}\nError: ${MONOFONIC_ERROR}")
endif()

# Check that output file was created
if(NOT EXISTS ${OUTPUT_FILE})
    message(FATAL_ERROR "monofonIC did not create output file: ${OUTPUT_FILE}")
endif()

message("monofonIC completed successfully")

# Run comparison script
message("Comparing output to reference: ${REFERENCE_FILE}")
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${COMPARE_SCRIPT} ${REFERENCE_FILE} ${OUTPUT_FILE}
    RESULT_VARIABLE COMPARE_RESULT
    OUTPUT_VARIABLE COMPARE_OUTPUT
    ERROR_VARIABLE COMPARE_ERROR
)

# Print comparison output
message("${COMPARE_OUTPUT}")
if(COMPARE_ERROR)
    message("${COMPARE_ERROR}")
endif()

# Check comparison result
if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR "Comparison failed: files differ")
endif()

message("Test passed: output matches reference")
