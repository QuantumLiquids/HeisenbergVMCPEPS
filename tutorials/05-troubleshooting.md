## Troubleshooting

This page is organized by failure stage:

1. Build-time issues
2. Runtime issues
3. Output/plot issues

### 1) Build Issues

### 1.1 OpenMP detection fails during CMake in macOS

Symptom:

- CMake fails to find OpenMP CXX flags/libs.

Reason:
 You use the apple clang rather llvm clang.


Action:
0. Home brew install the llvm clang.
1. Re-run CMake once with explicit llvm clang C and CXX compilers.
2. Use the same command again if first run still fails.

```bash
cmake -S . -B build \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
```

### 1.2 macOS libc++ ABI mismatch (`std::__1::__hash_memory`)

Symptom:

- Linker error like:
  - `Undefined symbols for architecture arm64: std::__1::__hash_memory`

Cause:

- Toolchain/header ABI mismatch between Homebrew LLVM and libc++ headers.

Action:

```bash
SDK="$(xcrun --sdk macosx --show-sdk-path)"
cmake -S . -B build \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
  -DCMAKE_OSX_SYSROOT="$SDK" \
  -DCMAKE_CXX_FLAGS="-nostdinc++ -isystem $SDK/usr/include/c++/v1"
```


### 2) Runtime Issues

### 2.1 Heavy `psi_consistency` warnings

Symptom:

- Many lines like `psi_consistency rel_err > threshold`.

Meaning:
- Please Read the paper for OBC PEPS

Typical mitigations:

- tighten contraction accuracy (`BMPSTruncErr` smaller, larger `Dbmps_max`, or TRG accuracy for PBC)

### 2.2 `sitps_tile` fails before dump

Common error classes:

- `--unit-ly/--unit-lx are not allowed for OBC mode`
- `PBC requires target dimensions to be integer multiples of unit cell`
- `source size >= 3x3` required for OBC boundary-preserving growth
- leg-consistency mismatch reports on virtual bonds

Actions:

1. Confirm source boundary condition from metadata:

```bash
cat tpsfinal/tps_meta.txt
```

2. For PBC, ensure target is a multiple of unit:

```bash
./sitps_tile \
  --input-dir tpsfinal \
  --output-dir tpsfinal_8x8 \
  --target-ly 8 \
  --target-lx 8 \
  --unit-ly 2 \
  --unit-lx 2
```

3. For OBC, remove unit flags and use source at least `3x3`.

4. If optional `--physics-json` is passed, ensure its `Lx/Ly/BoundaryCondition` exactly match the tiling target and SITPS BC.

### 2.3 PBC tiling succeeds but later VMC/TRG rejects size

Symptom:

- `sitps_tile` prints a warning about TRG compatibility.
- Later PBC VMC/measure aborts with TRG size constraints.

Reason:

- Current TRG backend expects square lattice and linear size in `L=2^k` or `L=3*2^k`.

Action:

- Choose target size satisfying TRG constraints for PBC runs.
- Or use OBC workflow if non-TRG sizes are required.

### 3) Output and Plot Issues

