# monofonIC Regression Tests

This directory contains regression tests for monofonIC to catch commits that break compatibility or change numerical results.

## Overview

The test suite consists of:
- **5 regression test configurations** covering different LPT orders, particle loads, and output formats
- **1 MPI consistency test** that verifies identical results across different MPI task counts
- **Reference HDF5 files** containing expected outputs
- **Comparison script** that performs hybrid tolerance checking (exact for integers, relative tolerance for floats)
- **CMake/CTest integration** for easy test execution
- **GitHub Actions CI** support for automated testing

## Test Cases

| Test Name | LPT Order | Particles | Baryons | Vrel | Output Format |
|-----------|-----------|-----------|---------|------|---------------|
| `test_1lpt_sc_generic` | 1LPT | sc (32³) | No | - | Generic HDF5 |
| `test_2lpt_sc_gadget` | 2LPT | sc (32³) | No | - | Gadget HDF5 |
| `test_3lpt_bcc_swift` | 3LPT | bcc (2×32³) | No | - | SWIFT HDF5 |
| `test_2lpt_baryons_generic` | 2LPT | sc (32³) | Yes | No | Generic HDF5 |
| `test_2lpt_baryons_vrel_gadget` | 2LPT | sc (32³) | Yes | Yes | Gadget HDF5 |

All tests use:
- Grid resolution: 32³ (fast execution)
- Box size: 100 Mpc/h
- Starting redshift: z = 50
- Transfer function: Eisenstein & Hu fitting formulae
- RNG: NGENIC with fixed seed (12345)
- Cosmology: Planck2018EE+BAO+SN

### MPI Consistency Test

| Test Name | Description |
|-----------|-------------|
| `test_mpi_consistency` | Runs the same 2LPT configuration with 1, 2, and 4 MPI tasks and verifies outputs are identical |

**Purpose**: Ensures that MPI parallelization is deterministic and doesn't introduce non-determinism or bugs. This test catches:
- MPI-related race conditions
- Domain decomposition errors
- Non-deterministic RNG behavior across MPI tasks
- Communication errors between MPI processes

**Requirements**:
- MPI must be enabled in the build (`ENABLE_MPI=ON`)
- `mpirun` or equivalent MPI launcher must be available

**Note**: This test is automatically skipped if MPI is not available.

## Running Tests

### Prerequisites

```bash
# Python 3 with h5py and numpy
pip3 install h5py numpy
```

### Building with Tests

```bash
# From repository root
mkdir build && cd build
cmake ..
make

# Generate reference files (only needed once, or after intentional changes)
bash ../tests/scripts/generate_references.sh

# Run all tests
ctest --output-on-failure

# Run specific test
ctest -R test_1lpt_sc_generic --verbose

# Run only MPI consistency test
ctest -R test_mpi_consistency --verbose

# Run only regression tests (exclude MPI test)
ctest -L regression -LE mpi --output-on-failure

# Run with parallel execution (note: doesn't speed up MPI test)
ctest -j4 --output-on-failure
```

### Manual Test Execution

You can also run tests manually:

```bash
cd build

# Run monofonIC with a test config
./monofonIC ../tests/configs/test_1lpt_sc_generic.conf

# Compare output to reference
python3 ../tests/scripts/compare_hdf5.py \
    ../tests/references/test_1lpt_sc_generic.hdf5 \
    test_1lpt_sc_generic.hdf5 \
    --verbose
```

## Comparison Logic

The `compare_hdf5.py` script uses **hybrid tolerance** for robust regression testing:

- **Integer datasets** (IDs, particle counts): Exact bit-for-bit comparison
- **Float datasets** (positions, velocities, masses): Relative tolerance (default: 1e-10)
- **Attributes**: Exact comparison

This approach ensures:
- Strict correctness for discrete quantities
- Robustness against minor floating-point variations across platforms/compilers
- Sensitive detection of actual physics/algorithm changes

### Tolerance Customization

```bash
# Use stricter tolerance
python3 scripts/compare_hdf5.py ref.hdf5 test.hdf5 --rtol 1e-12

# Use looser tolerance (not recommended)
python3 scripts/compare_hdf5.py ref.hdf5 test.hdf5 --rtol 1e-8
```

## Regenerating References

You may need to regenerate reference files when:
- Intentional algorithm improvements are made
- Physics corrections are implemented
- Output format changes are introduced

**Important**: Only regenerate references after carefully verifying that the changes are correct!

### Regeneration Process

```bash
cd build

# Option 1: Use the script (recommended)
bash ../tests/scripts/generate_references.sh

# Option 2: Manual regeneration
./monofonIC ../tests/configs/test_1lpt_sc_generic.conf
mv test_1lpt_sc_generic.hdf5 ../tests/references/

# Verify tests pass with new references
ctest --output-on-failure
```

### Committing New References

After regenerating references, commit them to the repository:

```bash
git add tests/references/*.hdf5
git commit -m "Update test references after [brief description of changes]"
```

Include in your commit message:
- Why references were updated
- What physical/numerical changes occurred
- Verification that results are correct

## Adding New Tests

To add a new regression test:

### 1. Create Configuration File

Create `tests/configs/test_mytest.conf`:

```ini
[setup]
GridRes         = 32
BoxLength       = 100
zstart          = 50.0
LPTorder        = 2
DoBaryons       = no
ParticleLoad    = sc

[cosmology]
ParameterSet    = Planck2018EE+BAO+SN
transfer        = eisenstein

[random]
generator       = NGENIC
seed            = 12345

[execution]
NumThreads      = 1

[output]
format          = generic
filename        = test_mytest.hdf5
generic_out_eulerian = no
```

### 2. Register Test in CMake

Edit `tests/CMakeLists.txt` and add:

```cmake
add_regression_test(
    test_mytest
    test_mytest.conf
    test_mytest.hdf5
)
```

### 3. Generate Reference

```bash
cd build
./monofonIC ../tests/configs/test_mytest.conf
mv test_mytest.hdf5 ../tests/references/
```

### 4. Verify Test

```bash
ctest -R test_mytest --verbose
```

### 5. Update Documentation

Update this README.md to document the new test case.

## Continuous Integration

Tests run automatically on GitHub Actions for every push and pull request to the `master` branch.

The CI workflow:
1. Installs dependencies (FFTW3, GSL, HDF5, Python3, h5py)
2. Builds monofonIC
3. Generates reference files
4. Runs all regression tests
5. Reports failures with detailed output

See `.github/workflows/cmake-multi-platform.yml` for details.

### CI Test Failure

If tests fail in CI:

1. Check the GitHub Actions logs for detailed error messages
2. The comparison script will show which datasets differ and by how much
3. Verify locally:
   ```bash
   git checkout <failing-commit>
   mkdir build && cd build
   cmake .. && make
   bash ../tests/scripts/generate_references.sh
   ctest --output-on-failure
   ```

## Test Design Philosophy

### Why Small Tests?

- **Fast execution**: 32³ grid completes in seconds
- **Frequent CI runs**: Developers get quick feedback
- **Low storage**: Reference files are small (~few MB total)
- **Still comprehensive**: Captures algorithm and output format correctness

### What Tests Catch

**Regression Tests:**
✓ Changes to particle positions/velocities
✓ Changes to LPT algorithm implementation
✓ Output format modifications
✓ Baryon physics changes
✓ RNG differences
✓ Cosmology calculation errors

**MPI Consistency Test:**
✓ MPI-related race conditions or non-determinism
✓ Domain decomposition errors
✓ MPI communication bugs
✓ Non-deterministic RNG across different MPI task counts

**Not Covered:**
✗ Performance regressions (not measured)
✗ Large-scale/convergence behavior (use full-scale tests)

## Troubleshooting

### Test fails with "Reference file not found"

```bash
cd build
bash ../tests/scripts/generate_references.sh
```

### Test fails with "h5py not found"

```bash
pip3 install h5py numpy
```

### Tests pass locally but fail in CI

- Check compiler differences (GCC version, flags)
- Verify dependency versions match
- Look for non-deterministic behavior (though tests use fixed seeds)

### Comparison shows small differences in floats

- If differences are at machine precision (~1e-15), this is expected
- Default tolerance (1e-10) should handle cross-platform variations
- If differences are larger, investigate the root cause before updating references

### MPI-related test failures

Current tests use `NumThreads = 1` and are designed for single-process execution.
For MPI testing, consider adding separate MPI-specific test cases.

## Future Enhancements

Potential improvements to the test suite:

- [ ] Add MPI-parallel test cases
- [ ] Test CLASS transfer function integration (separate from these fast tests)
- [ ] Add convergence tests (comparing different resolutions)
- [ ] Performance benchmarking framework
- [ ] Test PANPHASIA RNG (requires license agreement)
- [ ] Add glass initial conditions test
- [ ] Test mode fixing and phase inversion options

## Questions or Issues?

- Report test-related issues on GitHub: https://github.com/cosmo-sims/monofonIC/issues
- Join the MUSIC user group: https://groups.google.com/g/cosmo_music
