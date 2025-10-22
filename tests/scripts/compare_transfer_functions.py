#!/usr/bin/env python3
"""
Transfer function comparison script for monofonIC regression tests.

Compares two transfer function text files with relative tolerance.

Exit codes:
  0 - Files match within tolerance
  1 - Files differ
  2 - Error (file not found, cannot read, etc.)
"""

import sys
import argparse
import numpy as np


def load_transfer_file(filename):
    """
    Load transfer function data from a text file.

    Args:
        filename: Path to transfer function file

    Returns:
        dict: Dictionary with keys 'k' (array) and transfer function names (arrays)
    """
    try:
        # Read the file, skipping comment lines
        data = np.loadtxt(filename, comments='#')
    except Exception as e:
        print(f"ERROR: Could not read file '{filename}': {e}", file=sys.stderr)
        return None

    if data.shape[1] < 7:
        print(f"ERROR: File '{filename}' has insufficient columns (expected 7, got {data.shape[1]})",
              file=sys.stderr)
        return None

    # Parse columns
    result = {
        'k': data[:, 0],
        'delta_matter': data[:, 1],
        'delta_cdm': data[:, 2],
        'delta_baryon': data[:, 3],
        'theta_matter': data[:, 4],
        'theta_cdm': data[:, 5],
        'theta_baryon': data[:, 6],
    }

    return result


def compare_transfer_data(ref_data, test_data, rtol=1e-6, verbose=False):
    """
    Compare two transfer function datasets.

    Args:
        ref_data: Reference transfer function data (dict)
        test_data: Test transfer function data (dict)
        rtol: Relative tolerance for comparison
        verbose: Print detailed comparison results

    Returns:
        bool: True if all data matches within tolerance
    """
    # Check k values match
    if not np.allclose(ref_data['k'], test_data['k'], rtol=1e-10):
        print("ERROR: k values do not match between files", file=sys.stderr)
        max_k_diff = np.max(np.abs(ref_data['k'] - test_data['k']))
        print(f"       Max k difference: {max_k_diff:.3e}", file=sys.stderr)
        return False

    # Transfer function types to compare
    tf_types = ['delta_matter', 'delta_cdm', 'delta_baryon',
                'theta_matter', 'theta_cdm', 'theta_baryon']

    all_passed = True
    results = []

    for tf_type in tf_types:
        ref = ref_data[tf_type]
        test = test_data[tf_type]

        # Handle NaN values
        ref_nan = np.isnan(ref)
        test_nan = np.isnan(test)

        if not np.array_equal(ref_nan, test_nan):
            all_passed = False
            n_mismatch = np.sum(ref_nan != test_nan)
            results.append({
                'name': tf_type,
                'passed': False,
                'message': f"NaN mismatch ({n_mismatch} points)",
                'max_rel_diff': np.nan
            })
            continue

        # Compare non-NaN values
        valid_mask = ~ref_nan
        ref_valid = ref[valid_mask]
        test_valid = test[valid_mask]

        if len(ref_valid) == 0:
            # All values are NaN
            results.append({
                'name': tf_type,
                'passed': True,
                'message': "All values NaN",
                'max_rel_diff': 0.0
            })
            continue

        # Compute relative differences
        abs_ref = np.abs(ref_valid)
        abs_diff = np.abs(test_valid - ref_valid)

        # Use relative difference where |ref| > 1e-30, otherwise use absolute difference
        rel_diff = np.where(abs_ref > 1e-30, abs_diff / abs_ref, abs_diff)

        max_rel_diff = np.max(rel_diff)
        max_rel_diff_idx = np.argmax(rel_diff)

        # Find corresponding k value
        valid_indices = np.where(valid_mask)[0]
        k_at_max = ref_data['k'][valid_indices[max_rel_diff_idx]]

        # Check if within tolerance
        n_failures = np.sum(rel_diff > rtol)

        if n_failures > 0:
            all_passed = False
            results.append({
                'name': tf_type,
                'passed': False,
                'message': f"{n_failures}/{len(ref_valid)} points exceed tolerance",
                'max_rel_diff': max_rel_diff,
                'k_at_max': k_at_max
            })
        else:
            results.append({
                'name': tf_type,
                'passed': True,
                'message': "Passed",
                'max_rel_diff': max_rel_diff
            })

    # Print results
    if verbose or not all_passed:
        print()
        print("Transfer Function Comparison Results")
        print("=" * 80)
        print(f"Relative tolerance: {rtol:.3e}")
        print("-" * 80)

        for result in results:
            status = "PASS" if result['passed'] else "FAIL"
            print(f"{status:4s} | {result['name']:15s} | max_rel_diff = {result['max_rel_diff']:.3e} | {result['message']}")
            if not result['passed'] and 'k_at_max' in result:
                print(f"       Max difference at k = {result['k_at_max']:.3e} h/Mpc")

        print("-" * 80)

    if all_passed:
        print("SUCCESS: All transfer functions match within tolerance")
        return True
    else:
        n_failed = sum(1 for r in results if not r['passed'])
        print(f"FAILURE: {n_failed}/{len(results)} transfer functions failed comparison")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Compare two transfer function text files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s ref.txt test.txt
  %(prog)s ref.txt test.txt --rtol 1e-9
  %(prog)s ref.txt test.txt --verbose
        """
    )

    parser.add_argument('reference', help='Reference transfer function file')
    parser.add_argument('test', help='Test transfer function file')
    parser.add_argument('--rtol', type=float, default=1e-6,
                        help='Relative tolerance for comparison (default: 1e-6)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed comparison results')

    args = parser.parse_args()

    # Load reference file
    print(f"Loading reference: {args.reference}")
    ref_data = load_transfer_file(args.reference)
    if ref_data is None:
        return 2

    # Load test file
    print(f"Loading test file: {args.test}")
    test_data = load_transfer_file(args.test)
    if test_data is None:
        return 2

    # Compare
    if compare_transfer_data(ref_data, test_data, rtol=args.rtol, verbose=args.verbose):
        return 0
    else:
        return 1


if __name__ == '__main__':
    sys.exit(main())
