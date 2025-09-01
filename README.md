# finite-size Spin Models: Projected entangled pair state

___

## Models

- spin-$\frac{1}{2}$ square lattice $J_1$-$J_2$ Heisenberg and XY model
- spin-$\frac{1}{2}$ triangle lattice $J_1$-$J_2$ Heisenberg model

### Hamiltonian

The general $J_1$-$J_2$ Heisenberg model on a lattice is given by:
\[
H = J_1 \sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j + J_2 \sum_{\langle\langle i,j \rangle\rangle} \mathbf{S}_i \cdot \mathbf{S}_j
\]
where $\langle i,j \rangle$ denotes nearest-neighbor pairs and $\langle\langle i,j \rangle\rangle$ denotes next-nearest-neighbor pairs. $\mathbf{S}_i$ is the spin-$\frac{1}{2}$ operator at site $i$.

The general $J_1$-$J_2$ XY model on a lattice is given by:
\[
H = J_1 \sum_{\langle i,j \rangle} (S_i^x S_j^x+S_i^y S_j^y) + J_2 \sum_{\langle\langle i,j \rangle\rangle} (S_i^x S_j^x+S_i^y S_j^y).
\]

---

## Binaries and Their Functionality

After compiling (see Build), you obtain these binaries:

### Unified entries

| Binary | Function | Usage |
|---|---|---|
| simple_update | PEPS simple update (Square Heisenberg/XY, Triangle) | `./simple_update params/physics_params.json params/simple_update_algorithm_params.json` |
| vmc_optimize | VMC optimization | `mpirun -n 1 ./vmc_optimize params/physics_params.json params/vmc_algorithm_params.json` |
| mc_measure | Monte Carlo measurement | `./mc_measure params/physics_params.json params/measure_algorithm_params.json` |

### Notes

- Set `ModelType` in `physics_params.json` to dispatch: `SquareHeisenberg`, `SquareXY`, `TriangleHeisenberg`.
- Triangle uses square-lattice PEPS as ansatz under the hood.

---

## Dependence

- C++17+, CMake≥3.12
- BLAS/LAPACK: Intel MKL or OpenBLAS
- MPI
- Tensor: QuantumLiquids/TensorToolkit
- DMRG: QuantumLiquids/UltraDMRG
- PEPS: QuantumLiquids/PEPS

---

## Build Instructions

```bash
git clone https://github.com/QuantumLiquids/HeisenbergVMCPEPS.git
cd HeisenbergVMCPEPS
mkdir build && cd build
cmake ..
make -j4
```

Common options:
- `-DU1SYM=1` enable U(1) symmetry; default off
- `-DRealCode=ON` real tensors (default), OFF for complex

BLAS/LAPACK:
- x86: oneAPI MKL (single-thread) by default
- ARM: OpenBLAS by default

---

## Documentation

- [Tutorials](tutorials/) – short user guide and project plan
- API docs: see `qlten`/`qlpeps` upstream references (local docs/ if available)

## Quick Links

- [Quick Start](tutorials/01-quick-start.md)
- [Parameter System](tutorials/02-parameter-system.md)
- [Parameter Reference](tutorials/03-parameter-reference.md)
- [Workflows](tutorials/04-algorithm-workflows.md)
- [Optimizer Examples](tutorials/05-optimizer-examples.md)
 - Machine test: [square12x12D8vmc](tutorials/06-square12x12D8vmc.md)

---

## Author

Hao-Xin Wang — [wanghaoxin1996@gmail.com](mailto:wanghaoxin1996@gmail.com)
