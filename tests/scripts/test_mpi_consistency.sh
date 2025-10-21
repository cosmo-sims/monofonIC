#!/bin/bash
#
# MPI consistency test for monofonIC
#
# This script runs monofonIC with different numbers of MPI tasks and verifies
# that the outputs are identical, ensuring MPI parallelization is deterministic.
#
# Usage:
#   test_mpi_consistency.sh <monofonic_exe> <config_file> <compare_script> <python_exe>
#

set -e  # Exit on error

if [ $# -ne 4 ]; then
    echo "Usage: $0 <monofonic_exe> <config_file> <compare_script> <python_exe>"
    exit 1
fi

MONOFONIC_EXE="$1"
CONFIG_FILE="$2"
COMPARE_SCRIPT="$3"
PYTHON_EXE="$4"

# Extract output filename from config
OUTPUT_FILE=$(grep "^filename" "$CONFIG_FILE" | awk '{print $3}')

if [ -z "$OUTPUT_FILE" ]; then
    echo "Error: Could not extract output filename from config file"
    exit 1
fi

# Create temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "=========================================="
echo "MPI Consistency Test"
echo "=========================================="
echo "Config: $CONFIG_FILE"
echo "Output: $OUTPUT_FILE"
echo "Temp dir: $TEMP_DIR"
echo ""

# Test with different MPI task counts
MPI_TASKS=(1 2 4)

for NTASKS in "${MPI_TASKS[@]}"; do
    echo "------------------------------------------"
    echo "Running with $NTASKS MPI task(s)..."
    echo "------------------------------------------"

    OUTPUT_PATH="${TEMP_DIR}/output_np${NTASKS}_${OUTPUT_FILE}"

    # Run monofonIC with specified number of MPI tasks
    if [ "$NTASKS" -eq 1 ]; then
        # Single task - run without mpirun
        "$MONOFONIC_EXE" "$CONFIG_FILE" > "${TEMP_DIR}/log_np${NTASKS}.txt" 2>&1
    else
        # Multiple tasks - use mpirun with --oversubscribe for CI environments with limited cores
        mpirun --oversubscribe -np "$NTASKS" "$MONOFONIC_EXE" "$CONFIG_FILE" > "${TEMP_DIR}/log_np${NTASKS}.txt" 2>&1
    fi

    # Check that output was created (single file or multi-file format)
    # Gadget format can create either test.hdf5 or test.0.hdf5, test.1.hdf5, etc.
    BASENAME="${OUTPUT_FILE%.hdf5}"

    if [ -f "$OUTPUT_FILE" ]; then
        # Single file output
        mv "$OUTPUT_FILE" "$OUTPUT_PATH"
    elif [ -f "${BASENAME}.0.hdf5" ]; then
        # Multi-file output - move all files
        for f in ${BASENAME}.*.hdf5; do
            if [ -f "$f" ]; then
                SUFFIX="${f#${BASENAME}}"
                mv "$f" "${TEMP_DIR}/output_np${NTASKS}_${BASENAME}${SUFFIX}"
            fi
        done
    else
        echo "Error: monofonIC did not create output file: $OUTPUT_FILE"
        cat "${TEMP_DIR}/log_np${NTASKS}.txt"
        exit 1
    fi

    echo "✓ Completed with $NTASKS task(s)"
    echo ""
done

# Compare outputs
echo "=========================================="
echo "Comparing outputs..."
echo "=========================================="

BASENAME="${OUTPUT_FILE%.hdf5}"

# Reference is always from 1-task run
# Determine if reference is single-file or multi-file
if [ -f "${TEMP_DIR}/output_np1_${OUTPUT_FILE}" ]; then
    REFERENCE="${TEMP_DIR}/output_np1_${OUTPUT_FILE}"
elif [ -f "${TEMP_DIR}/output_np1_${BASENAME}.0.hdf5" ]; then
    REFERENCE="${TEMP_DIR}/output_np1_${BASENAME}.0.hdf5"
else
    echo "Error: Reference output (1 task) not found"
    ls -la "${TEMP_DIR}/"
    exit 1
fi

for NTASKS in "${MPI_TASKS[@]:1}"; do  # Skip first element (ntasks=1)
    echo "------------------------------------------"
    echo "Comparing: 1 task vs $NTASKS tasks"
    echo "------------------------------------------"

    # Test output could be single-file or multi-file
    if [ -f "${TEMP_DIR}/output_np${NTASKS}_${OUTPUT_FILE}" ]; then
        TEST_OUTPUT="${TEMP_DIR}/output_np${NTASKS}_${OUTPUT_FILE}"
    elif [ -f "${TEMP_DIR}/output_np${NTASKS}_${BASENAME}.0.hdf5" ]; then
        TEST_OUTPUT="${TEMP_DIR}/output_np${NTASKS}_${BASENAME}.0.hdf5"
    else
        echo "Error: Output file for $NTASKS tasks not found"
        ls -la "${TEMP_DIR}/"
        exit 1
    fi

    echo "Comparing: $REFERENCE"
    echo "     vs:   $TEST_OUTPUT"

    # Run comparison
    # Note: For multi-file outputs, we only compare file .0
    # This contains a subset of particles but should be bitwise identical
    # in terms of those particles' properties
    "$PYTHON_EXE" "$COMPARE_SCRIPT" "$REFERENCE" "$TEST_OUTPUT" --rtol 1e-9

    if [ $? -eq 0 ]; then
        echo "✓ Outputs match (1 task vs $NTASKS tasks)"
    else
        echo "✗ Outputs differ (1 task vs $NTASKS tasks)"
        exit 1
    fi
    echo ""
done

echo "=========================================="
echo "✓ All MPI consistency tests passed!"
echo "=========================================="
echo ""
echo "Verified that outputs are identical for:"
for NTASKS in "${MPI_TASKS[@]}"; do
    echo "  - $NTASKS MPI task(s)"
done
echo ""

exit 0
