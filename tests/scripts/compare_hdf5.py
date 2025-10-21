#!/usr/bin/env python3
"""
HDF5 comparison script for monofonIC regression tests.

Compares two HDF5 files with hybrid tolerance:
- Integer datasets: exact bit-for-bit comparison
- Float datasets: relative tolerance comparison (default 1e-10)

Exit codes:
  0 - Files match within tolerance
  1 - Files differ
  2 - Error (file not found, cannot read, etc.)
"""

import sys
import argparse
import numpy as np
import h5py


def is_integer_dtype(dtype):
    """Check if numpy dtype is an integer type."""
    return np.issubdtype(dtype, np.integer)


def is_float_dtype(dtype):
    """Check if numpy dtype is a floating point type."""
    return np.issubdtype(dtype, np.floating)


def compare_datasets(dataset1, dataset2, name, rtol=1e-10, atol=1e-14):
    """
    Compare two HDF5 datasets with appropriate tolerance.

    Args:
        dataset1: First HDF5 dataset
        dataset2: Second HDF5 dataset
        name: Name of the dataset (for error reporting)
        rtol: Relative tolerance for float comparison
        atol: Absolute tolerance for float comparison

    Returns:
        (bool, str): (True if match, error message if no match)
    """
    # Check shapes match
    if dataset1.shape != dataset2.shape:
        return False, f"Shape mismatch: {dataset1.shape} vs {dataset2.shape}"

    # Check dtypes match
    if dataset1.dtype != dataset2.dtype:
        return False, f"Dtype mismatch: {dataset1.dtype} vs {dataset2.dtype}"

    # Read data
    data1 = dataset1[...]
    data2 = dataset2[...]

    # Integer comparison: exact match
    if is_integer_dtype(dataset1.dtype):
        if not np.array_equal(data1, data2):
            diff_mask = data1 != data2
            num_diff = np.sum(diff_mask)
            return False, f"Integer data mismatch: {num_diff}/{data1.size} elements differ"
        return True, ""

    # Float comparison: relative tolerance
    elif is_float_dtype(dataset1.dtype):
        if not np.allclose(data1, data2, rtol=rtol, atol=atol):
            # Compute relative and absolute differences
            abs_diff = np.abs(data1 - data2)
            rel_diff = abs_diff / (np.abs(data1) + atol)
            max_abs_diff = np.max(abs_diff)
            max_rel_diff = np.max(rel_diff)
            return False, (f"Float data mismatch: max_abs_diff={max_abs_diff:.3e}, "
                          f"max_rel_diff={max_rel_diff:.3e} (rtol={rtol:.3e})")
        return True, ""

    # String or other types: exact comparison
    else:
        if not np.array_equal(data1, data2):
            return False, "Data mismatch (non-numeric)"
        return True, ""


def compare_attributes(attrs1, attrs2, name):
    """
    Compare attributes of two HDF5 objects.

    Args:
        attrs1: First attribute dict
        attrs2: Second attribute dict
        name: Name of the object (for error reporting)

    Returns:
        (bool, str): (True if match, error message if no match)
    """
    # Check all keys present in both
    keys1 = set(attrs1.keys())
    keys2 = set(attrs2.keys())

    if keys1 != keys2:
        missing_in_2 = keys1 - keys2
        missing_in_1 = keys2 - keys1
        msg = ""
        if missing_in_2:
            msg += f"Missing in file2: {missing_in_2}. "
        if missing_in_1:
            msg += f"Missing in file1: {missing_in_1}."
        return False, msg

    # Compare each attribute
    for key in keys1:
        val1 = attrs1[key]
        val2 = attrs2[key]

        # Handle scalar vs array
        if np.isscalar(val1):
            val1 = np.array([val1])
        if np.isscalar(val2):
            val2 = np.array([val2])

        # Convert to arrays for comparison
        val1 = np.asarray(val1)
        val2 = np.asarray(val2)

        if val1.dtype != val2.dtype:
            return False, f"Attribute '{key}' dtype mismatch: {val1.dtype} vs {val2.dtype}"

        if not np.array_equal(val1, val2):
            return False, f"Attribute '{key}' value mismatch: {val1} vs {val2}"

    return True, ""


def compare_groups_recursive(group1, group2, path="/", rtol=1e-10, atol=1e-14, verbose=False):
    """
    Recursively compare two HDF5 groups.

    Args:
        group1: First HDF5 group
        group2: Second HDF5 group
        path: Current path in HDF5 hierarchy
        rtol: Relative tolerance for float comparison
        atol: Absolute tolerance for float comparison
        verbose: Print detailed comparison info

    Returns:
        (bool, list): (True if match, list of error messages)
    """
    errors = []

    # Compare attributes
    match, msg = compare_attributes(group1.attrs, group2.attrs, path)
    if not match:
        errors.append(f"{path} [attributes]: {msg}")

    # Check all keys present in both groups
    keys1 = set(group1.keys())
    keys2 = set(group2.keys())

    if keys1 != keys2:
        missing_in_2 = keys1 - keys2
        missing_in_1 = keys2 - keys1
        if missing_in_2:
            errors.append(f"{path}: Missing in file2: {missing_in_2}")
        if missing_in_1:
            errors.append(f"{path}: Missing in file1: {missing_in_1}")
        # Continue with common keys
        keys_common = keys1 & keys2
    else:
        keys_common = keys1

    # Compare each item
    for key in sorted(keys_common):
        item1 = group1[key]
        item2 = group2[key]
        item_path = f"{path}{key}" if path == "/" else f"{path}/{key}"

        # Check if both are the same type (group or dataset)
        if isinstance(item1, h5py.Group) and isinstance(item2, h5py.Group):
            if verbose:
                print(f"Comparing group: {item_path}")
            _, sub_errors = compare_groups_recursive(item1, item2, item_path + "/",
                                                     rtol, atol, verbose)
            errors.extend(sub_errors)

        elif isinstance(item1, h5py.Dataset) and isinstance(item2, h5py.Dataset):
            if verbose:
                print(f"Comparing dataset: {item_path}")
            match, msg = compare_datasets(item1, item2, item_path, rtol, atol)
            if not match:
                errors.append(f"{item_path}: {msg}")

            # Also compare dataset attributes
            match, msg = compare_attributes(item1.attrs, item2.attrs, item_path)
            if not match:
                errors.append(f"{item_path} [attributes]: {msg}")

        else:
            errors.append(f"{item_path}: Type mismatch (group vs dataset)")

    return len(errors) == 0, errors


def compare_hdf5_files(file1_path, file2_path, rtol=1e-10, atol=1e-14, verbose=False):
    """
    Compare two HDF5 files.

    Args:
        file1_path: Path to first HDF5 file
        file2_path: Path to second HDF5 file
        rtol: Relative tolerance for float comparison
        atol: Absolute tolerance for float comparison
        verbose: Print detailed comparison info

    Returns:
        (bool, list): (True if files match, list of error messages)
    """
    try:
        with h5py.File(file1_path, 'r') as f1, h5py.File(file2_path, 'r') as f2:
            return compare_groups_recursive(f1, f2, rtol=rtol, atol=atol, verbose=verbose)
    except FileNotFoundError as e:
        return False, [f"File not found: {e}"]
    except Exception as e:
        return False, [f"Error reading HDF5 files: {e}"]


def main():
    parser = argparse.ArgumentParser(
        description="Compare two HDF5 files for monofonIC regression testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Comparison strategy:
  - Integer datasets: exact bit-for-bit comparison
  - Float datasets: relative tolerance comparison (default rtol=1e-10)
  - Attributes: exact comparison

Exit codes:
  0 - Files match within tolerance
  1 - Files differ
  2 - Error (file not found, cannot read, etc.)
        """
    )

    parser.add_argument("file1", help="First HDF5 file (reference)")
    parser.add_argument("file2", help="Second HDF5 file (test output)")
    parser.add_argument("--rtol", type=float, default=1e-10,
                       help="Relative tolerance for float comparison (default: 1e-10)")
    parser.add_argument("--atol", type=float, default=1e-14,
                       help="Absolute tolerance for float comparison (default: 1e-14)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Print detailed comparison information")

    args = parser.parse_args()

    if args.verbose:
        print(f"Comparing HDF5 files:")
        print(f"  Reference: {args.file1}")
        print(f"  Test:      {args.file2}")
        print(f"  Tolerances: rtol={args.rtol:.3e}, atol={args.atol:.3e}")
        print()

    match, errors = compare_hdf5_files(args.file1, args.file2,
                                       rtol=args.rtol, atol=args.atol,
                                       verbose=args.verbose)

    if match:
        print(f"✓ Files match within tolerance (rtol={args.rtol:.3e})")
        return 0
    else:
        print(f"✗ Files differ! Found {len(errors)} difference(s):")
        for i, error in enumerate(errors, 1):
            print(f"  {i}. {error}")

        # Check if it's a real error (file not found, etc.)
        if any("not found" in err.lower() or "error reading" in err.lower()
               for err in errors):
            return 2
        return 1


if __name__ == "__main__":
    sys.exit(main())
