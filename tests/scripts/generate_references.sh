#!/bin/bash
#
# Script to generate reference HDF5 files for monofonIC regression tests
#
# Usage:
#   From build directory: cd build && bash ../tests/scripts/generate_references.sh
#   Or from tests directory in build: cd build/tests && bash ../../tests/scripts/generate_references.sh
#

set -e  # Exit on error

# Determine script directory (source tree)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_SOURCE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="${TESTS_SOURCE_DIR}/configs"
REFERENCE_DIR="${TESTS_SOURCE_DIR}/references"

# Find monofonIC executable
# Try different possible locations
MONOFONIC_EXE=""
if [ -f "./monofonIC" ]; then
    MONOFONIC_EXE="./monofonIC"
elif [ -f "../monofonIC" ]; then
    MONOFONIC_EXE="../monofonIC"
elif [ -f "../../monofonIC" ]; then
    MONOFONIC_EXE="../../monofonIC"
else
    echo "Error: Cannot find monofonIC executable"
    echo "Please run this script from the build directory or build/tests directory"
    echo "Example: cd build && bash ../tests/scripts/generate_references.sh"
    exit 1
fi

MONOFONIC_EXE="$(cd "$(dirname "$MONOFONIC_EXE")" && pwd)/$(basename "$MONOFONIC_EXE")"
echo "Using monofonIC executable: $MONOFONIC_EXE"

# Create references directory if it doesn't exist
mkdir -p "$REFERENCE_DIR"

# List of test configurations
TESTS=(
    "test_1lpt_sc_generic.conf:test_1lpt_sc_generic.hdf5"
    "test_2lpt_sc_gadget.conf:test_2lpt_sc_gadget.hdf5"
    "test_3lpt_bcc_swift.conf:test_3lpt_bcc_swift.hdf5"
    "test_2lpt_baryons_generic.conf:test_2lpt_baryons_generic.hdf5"
    "test_2lpt_baryons_vrel_gadget.conf:test_2lpt_baryons_vrel_gadget.hdf5"
)

# Temporary directory for running tests
TEMP_DIR=$(mktemp -d)
cd "$TEMP_DIR"

echo ""
echo "Generating reference outputs in temporary directory: $TEMP_DIR"
echo "Reference files will be saved to: $REFERENCE_DIR"
echo ""

# Generate each reference
for TEST_SPEC in "${TESTS[@]}"; do
    # Split on colon
    IFS=':' read -r CONFIG_FILE OUTPUT_FILE <<< "$TEST_SPEC"

    CONFIG_PATH="${CONFIG_DIR}/${CONFIG_FILE}"

    echo "=========================================="
    echo "Generating: $OUTPUT_FILE"
    echo "Config: $CONFIG_FILE"
    echo "=========================================="

    # Check config exists
    if [ ! -f "$CONFIG_PATH" ]; then
        echo "Error: Config file not found: $CONFIG_PATH"
        exit 1
    fi

    # Run monofonIC
    "$MONOFONIC_EXE" "$CONFIG_PATH"

    # Check output was created
    if [ ! -f "$OUTPUT_FILE" ]; then
        echo "Error: monofonIC did not create expected output: $OUTPUT_FILE"
        exit 1
    fi

    # Move to reference directory
    mv "$OUTPUT_FILE" "$REFERENCE_DIR/"
    echo "✓ Saved: ${REFERENCE_DIR}/${OUTPUT_FILE}"
    echo ""
done

# Cleanup
cd -
rm -rf "$TEMP_DIR"

echo "=========================================="
echo "All reference files generated successfully!"
echo "=========================================="
echo ""
echo "Reference files are now in: $REFERENCE_DIR"
echo ""
echo "You can now run the regression tests with:"
echo "  cd build"
echo "  ctest --output-on-failure"
echo ""
echo "Or run specific tests:"
echo "  ctest -R test_1lpt_sc_generic --verbose"
echo ""
