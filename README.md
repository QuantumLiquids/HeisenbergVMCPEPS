# finite-size Spin Models: Projected entangled pair state


___

## Models

- spin-$\frac{1}{2}$ square lattice $J_1-J_2$ Heisenberg model
- spin-$\frac{1}{2}$ triangle lattice $J_1-J_2$ Heisenberg model

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

1. Clone the repository:
   ```bash
   git clone https://github.com/QuantumLiquids/HeisenbergVMCPEPS.git
   cd HeisenbergVMCPEPS
   ```

2. Build using CMake:
   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_CXX_COMPILER=your_cxx_compiler \
            -DCMAKE_PREFIX_PATH=hint_to_find_tensortoolkit_ultradmrg_and_peps \
            -DU1SYM=1_if_want_use_tensor_u1_sym

   make -j4 
   ```

### BLAS/LAPACK Setting:
- On **x86 platforms**, CMakeList.txt will automatically use oneAPI single-threaded MKL as the BLAS implementation.
- On **ARM platforms**, it will default to OpenBLAS.
- If you want to link against a different implementation, such as an older version of MKL, sequential MKL, or AOCL, you will need to manually set the following CMake options:
    - *BLA_VENDOR*: [CMake Manual](https://cmake.org/cmake/help/latest/module/FindBLAS.html)
    - *BLAS_INCLUDE_DIR*: Set the header paths for the chosen BLAS/LAPACK implementation.

---

## Author

For inquiries, questions, or collaboration opportunities, please contact Hao-Xin via email:
[wanghaoxin1996@gmail.com](mailto:wanghaoxin1996@gmail.com).

## Acknowledgments
