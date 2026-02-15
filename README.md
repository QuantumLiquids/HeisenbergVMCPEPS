# finite-size Spin Models: Projected Entangled Pair States

This repository provides user-facing drivers for finite-size PEPS workflows on top of `qlpeps`.

## Models

- Spin-1/2 square lattice `J1-J2` Heisenberg model
- Spin-1/2 square lattice `J1-J2` XY model
- Spin-1/2 triangle lattice `J1-J2` Heisenberg model

## Binaries and Scope

Main workflow binaries:

| Binary | Purpose | Typical command |
|---|---|---|
| `simple_update` | Prepare initial PEPS/SITPS state | `./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json` |
| `vmc_optimize` | Optimize SITPS with VMC | `mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json` |
| `mc_measure` | Run Monte Carlo measurement on SITPS | `mpirun -n 1 ./mc_measure ../params/physics_params.json ../params/measure_algorithm_params.json` |
| `sitps_tile` | Tile/replicate SITPS to larger lattices with consistency checks | `./sitps_tile --input-dir tpsfinal --output-dir tpsfinal_8x8 --target-ly 8 --target-lx 8 --unit-ly 2 --unit-lx 2` |

DMRG benchmark binaries (appendix scope):

- `triangle_vmps`
- `triangle_dmrg`
- `triangle_dmrg_measure`
- `triangle_dmrg_measure_sf`

Notes:

- `ModelType` controls solver dispatch: `SquareHeisenberg`, `SquareXY`, `TriangleHeisenberg`.
- Triangle PBC is not supported in current PEPS backend.
- `src/kagome*` code is deprecated and corresponding binaries are disabled in default CMake.

## Dependencies

- C++20 compiler
- CMake >= 3.14
- BLAS/LAPACK (Intel MKL or OpenBLAS)
- MPI
- OpenMP (CXX)
- QuantumLiquids/TensorToolkit
- QuantumLiquids/UltraDMRG
- QuantumLiquids/PEPS

## Build

Generic:

```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
mkdir -p build
cmake -S . -B build
cmake --build build -j4
```

macOS/Homebrew LLVM (recommended when you hit libc++ ABI mismatch):

```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
./scripts/configure_macos_llvm_sdk.sh build_llvm_sdk Release
cmake --build build_llvm_sdk -j4
```

## Documentation (Read in Order)

Start here for the full user workflow:

1. `tutorials/01-quick-start.md` (run this)
2. `tutorials/02-concepts.md` (understand this)
3. `tutorials/03-recipes.md` (copy/tune this)
4. `tutorials/04-parameter-reference.md` (reference this)
5. `tutorials/05-troubleshooting.md` (debug this)
6. `tutorials/06-plotting-workflow.md` (plot this)
7. `tutorials/appendix-dmrg-benchmark.md` (benchmark appendix)

Legacy snapshots remain under `tutorials/deprecated/`.

## Author

Hao-Xin Wang — [wanghaoxin1996@gmail.com](mailto:wanghaoxin1996@gmail.com)
