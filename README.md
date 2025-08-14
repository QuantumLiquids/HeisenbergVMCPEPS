# finite-size Spin Models: Projected entangled pair state

___

## Models

- spin-$\frac{1}{2}$ square lattice $J_1$-$J_2$ Heisenberg and XY model
- spin-$\frac{1}{2}$ triangle lattice $J_1$-$J_2$ Heisenberg model

### Hamiltonian

The general $J_1$-$J_2$ Heisenberg model on a lattice is given by:
$$
H = J_1 \sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j + J_2 \sum_{\langle\langle i,j \rangle\rangle} \mathbf{S}_i \cdot \mathbf{S}_j
$$
where $\langle i,j \rangle$ denotes nearest-neighbor pairs and $\langle\langle i,j \rangle\rangle$ denotes
next-nearest-neighbor pairs. $\mathbf{S}_i$ is the spin-$\frac{1}{2}$ operator at site $i$.

The general $J_1$-$J_2$ XY model on a lattice is given by:
$$
H = J_1 \sum_{\langle i,j \rangle} ({S}_i^x \cdot {S}_j^x+{S}_i^y \cdot {S}_j^y) + J_2 \sum_{\langle\langle i,j \rangle\rangle} ({S}_i^x \cdot {S}_j^x+{S}_i^y \cdot {S}_j^y).
$$


---

## Binaries and Their Functionality
After compiling the programs according the **Build Instructions** below, 
you can obtain the following binaries with specific functionality.
### $J_1$-$J_2$ Heisenberg Model on Square Lattice

| Binary Name              | Functionality                                   | Usage                                                                        |
|--------------------------|-------------------------------------------------|------------------------------------------------------------------------------|
| **square_simple_update** | Simple update algorithm for PEPS.               | ./square_simple_update simple_update_params.json                             |
| **square_vmc_update**    | Variational Monte Carlo (VMC) update for  PEPS. | mpirun -np ${num_of_markov_chains} ./square_vmc_update vmc_update_params.json |
| **square_mc_measure**    | Monte Carlo measurement routines for PEPS.      | mpirun -np ${num_of_markov_chains} ./square_mc_measure vmc_update_params.json |

### $J_1$-$J_2$ XY Model on Square Lattice

| Binary Name                     | Functionality           | Usage                                                                                |
|---------------------------------|-------------------------|--------------------------------------------------------------------------------------|
| **square_planar_simple_update** | Simple update for PEPS. | ./square_planar_simple_update simple_update_params.json                              |
| **square_planar_vmc_update**    | VMC update for PEPS.    | mpirun -np ${num_of_markov_chains} ./square_planar_vmc_update vmc_update_params.json |

### $J_1$-$J_2$ Heisenberg Model on Triangle Lattice

For the triangle model calculation by PEPS, we use square-lattice PEPS as ansatz.

| Binary Name                   | Functionality                                                                      | Usage          |
|-------------------------------|------------------------------------------------------------------------------------|----------------|
| **triangle_simple_update**    | Simple update for PEPS.                                                            | Same as square |
| **triangle_loop_update**      | Loop update for  PEPS.                                                             | Same as square |
| **triangle_vmc_update**       | VMC update for PEPS (with total$S^z =0$ conservation in the Monte-Carlo sampling). | Same as square |
| **triangle_vmc_update_no_u1** | VMC update for PEPS (without $S^z$ conservation in the Monte-Carlo sampling).      | Same as square |
| **triangle_mc_measure**       | Monte Carlo measurement for  PEPS.                                                 | Same as square |
| **triangle_vmps**             | MPO-based DMRG for the model.                                                      |                |
| **triangle_dmrg**             | Traditional DMRG for the model.                                                    |                |
| **triangle_dmrg_measure**     | DMRG measurement routines for the model.                                           |                |
| **triangle_dmrg_measure_sf**  | DMRG measurement for spin fluctuation for the model.                               |                |

### Other Utilities

| Binary Name     | Functionality                                                                        |
|-----------------|--------------------------------------------------------------------------------------|
| **increaseD**   | Increase bond dimension $D$ for tensor network states.                               |
| **extend**      | System size extension for TPS.                                                       |
| **zero_evolve** | Imaginary time evolution for square lattice PEPS with $\tau=0$ to canonicalize PEPS. |

---

## Machine Test

### Square 12x12 D8 VMC Test

To run the machine test located in `machine_test/square12x12D8vmc/`, use code compiled **without U1 symmetric**. This means:

 Set `-DU1SYM=0` or leave the U1SYM option unset during CMake configuration

---

## Dependence

- **C++ Compiler:** C++17 or above
- **Build System:** CMake (version 3.12 or higher)
- **Math Libraries:** Intel MKL or OpenBLAS
- **Parallelization:** MPI
- **Tensor:** [QuantumLiquids/TensorToolkit](https://github.com/QuantumLiquids/TensorToolkit)
- **DMRG:** [QuantumLiquids/UltraDMRG](https://github.com/QuantumLiquids/UltraDMRG)
- **PEPS:** [QuantumLiquids/PEPS](https://github.com/QuantumLiquids/PEPS)

___

## Build Instructions

To compile the project using CMake, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/QuantumLiquids/HeisenbergVMCPEPS.git
   cd HeisenbergVMCPEPS
   ```

2. **Configure and build:**
   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_CXX_COMPILER=your_cxx_compiler \
            -DCMAKE_PREFIX_PATH=path_to_dependencies(tensor, dmrg, peps) \
            [additional options]
   make -j4
   ```

   **Common CMake options:**

    - `-DU1SYM=1` : Use U(1) symmetry-conserving tensor code  
      `-DU1SYM=0` or unset: No symmetry constraint (default)
    - `-DRealCode=ON` : Use real number tensors (default)  
      `-DRealCode=OFF` : Use complex number tensors

   **BLAS/LAPACK settings:**

    - On **x86 platforms**, CMakeLists.txt will automatically use oneAPI single-threaded MKL as the BLAS implementation.
    - On **ARM platforms**, it will default to OpenBLAS.
    - To use a different BLAS/LAPACK implementation, set:
        - `-DBLA_VENDOR=YourVendor` ([CMake Manual](https://cmake.org/cmake/help/latest/module/FindBLAS.html))
        - `-DBLAS_INCLUDE_DIR=path_to_blas_headers`

---

## Author

For inquiries, questions, or collaboration opportunities, please contact Hao-Xin via email:
[wanghaoxin1996@gmail.com](mailto:wanghaoxin1996@gmail.com)
