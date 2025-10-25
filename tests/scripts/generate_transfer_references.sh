#!/bin/bash
#
# Script to generate reference text files for monofonIC transfer function regression tests
#
# Usage:
#   From build directory: cd build && bash ../tests/scripts/generate_transfer_references.sh
#   Or from tests directory in build: cd build/tests && bash ../../tests/scripts/generate_transfer_references.sh
#

set -e  # Exit on error

# Determine script directory (source tree)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_SOURCE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="${TESTS_SOURCE_DIR}/configs/transfer"
REFERENCE_DIR="${TESTS_SOURCE_DIR}/references/transfer"

# Find test_transfer_functions executable
# Try different possible locations
TEST_EXE=""
if [ -f "./test_transfer_functions" ]; then
    TEST_EXE="./test_transfer_functions"
elif [ -f "./tests/test_transfer_functions" ]; then
    TEST_EXE="./tests/test_transfer_functions"
elif [ -f "../test_transfer_functions" ]; then
    TEST_EXE="../test_transfer_functions"
else
    echo "Error: Cannot find test_transfer_functions executable"
    echo "Please run this script from the build directory or build/tests directory"
    echo "Example: cd build && bash ../tests/scripts/generate_transfer_references.sh"
    exit 1
fi

TEST_EXE="$(cd "$(dirname "$TEST_EXE")" && pwd)/$(basename "$TEST_EXE")"
echo "Using test executable: $TEST_EXE"

# Create references directory structure if it doesn't exist
mkdir -p "$REFERENCE_DIR"

# List of plugins to test
PLUGINS=(
    "CLASS"
    "eisenstein"
)

# List of test configurations
# Format: config_name:output_basename
CONFIGS=(
    "fiducial.conf:fiducial.txt"
    "dark_energy.conf:dark_energy.txt"
    "massive_nu.conf:massive_nu.txt"
    "low_omega_m.conf:low_omega_m.txt"
    "high_z.conf:high_z.txt"
)

# Function to check if a config should be skipped for a plugin
should_skip_config() {
    local plugin=$1
    local config=$2

    # No special skipping needed currently
    return 1  # false, should not skip
}

echo ""
echo "Generating transfer function reference files"
echo "Reference files will be saved to: $REFERENCE_DIR"
echo ""

# Generate references for each plugin and configuration
for PLUGIN in "${PLUGINS[@]}"; do
    echo "========================================"
    echo "Plugin: $PLUGIN"
    echo "========================================"

    # Create plugin-specific directory
    PLUGIN_REF_DIR="${REFERENCE_DIR}/${PLUGIN}"
    mkdir -p "$PLUGIN_REF_DIR"

    for CONFIG_SPEC in "${CONFIGS[@]}"; do
        # Split on colon
        IFS=':' read -r CONFIG_FILE OUTPUT_FILE <<< "$CONFIG_SPEC"

        # Check if this config should be skipped for this plugin
        if should_skip_config "$PLUGIN" "$CONFIG_FILE"; then
            echo "  Skipping $CONFIG_FILE (not supported by $PLUGIN)"
            continue
        fi

        CONFIG_PATH="${CONFIG_DIR}/${CONFIG_FILE}"
        OUTPUT_PATH="${PLUGIN_REF_DIR}/${OUTPUT_FILE}"

        echo "  Generating: ${PLUGIN}/${OUTPUT_FILE}"
        echo "    Config: $CONFIG_FILE"

        # Check config exists
        if [ ! -f "$CONFIG_PATH" ]; then
            echo "    Error: Config file not found: $CONFIG_PATH"
            exit 1
        fi

        # Run test program in generate mode
        "$TEST_EXE" --generate "$PLUGIN" "$CONFIG_PATH" "$OUTPUT_PATH"

        # Check output was created
        if [ ! -f "$OUTPUT_PATH" ]; then
            echo "    Error: Test program did not create expected output: $OUTPUT_PATH"
            exit 1
        fi

        echo "    ✓ Saved: $OUTPUT_PATH"
    done

    echo ""
done

echo "=========================================="
echo "All transfer function reference files generated successfully!"
echo "=========================================="
echo ""
echo "Reference directory structure:"
find "$REFERENCE_DIR" -type f -name "*.txt" | sort
echo ""
echo "You can now run the transfer function regression tests with:"
echo "  cd build"
echo "  ctest -R test_transfer --output-on-failure"
echo ""
echo "Or run specific tests:"
echo "  ctest -R test_transfer_CLASS_fiducial --verbose"
echo ""
