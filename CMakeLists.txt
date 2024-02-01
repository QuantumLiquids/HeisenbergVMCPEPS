cmake_minimum_required(VERSION 3.14)
project(HeisenbergVMCPEPS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/")

# Build type.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

option(QLTEN_TIMING_MODE "Timing mode." OFF)
if (QLTEN_TIMING_MODE)
    add_definitions(-DQLTEN_TIMING_MODE)
endif ()

option(QLMPS_TIMING_MODE "Timing mode." OFF)
if (QLMPS_TIMING_MODE)
    add_definitions(-DQLMPS_TIMING_MODE)
endif ()

option(QLPEPS_TIMING_MODE "Timing mode." OFF)
if (QLPEPS_TIMING_MODE)
    add_definitions(-DQLPEPS_TIMING_MODE)
endif ()

option(QLTEN_MPI_TIMING_MODE "Timing mode." OFF)
if (QLTEN_MPI_TIMING_MODE)
    add_definitions(-DQLTEN_MPI_TIMING_MODE)
endif ()

option(QLPEPS_MPI_TIMING_MODE "Timing mode." OFF)
if (QLPEPS_MPI_TIMING_MODE)
    add_definitions(-DQLPEPS_MPI_TIMING_MODE)
endif ()

option(QLMPS_MPI_TIMING_MODE "Timing mode." OFF)
if (QLMPS_MPI_TIMING_MODE)
    add_definitions(-DQLMPS_MPI_TIMING_MODE)
endif ()

option(COMPILE_FIX_CODE "compile the fixing code" OFF)

# Compilation and linking control.
add_definitions(-Wall -g)

if (U1SYM)
    add_definitions(-DU1SYM)
    message("USE U1 symmetry code.")
endif ()

if (APPLE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        #Need test
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # May not work
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Note as of Intel oneAPI 2021.2, "source /opt/intel/oneapi/compiler/latest/env/vars.sh"
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -Wl, -rpath $ENV{CMPLR_ROOT}/mac/compiler/lib/libiomp5.dylib -liomp5 -lpthread -lm -ldl)
    endif ()
elseif (UNIX)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/.local/")
    message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

    #UNIX include APPLE, but we except it here
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Link the Intel's OpenMP library to avoid performance issue when the library calls the MKL's gesdd function.
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -Wl,--no-as-needed -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
    endif ()
endif ()


find_path(TENSOR_HEADER_PATH "qlten")
find_path(MPS_HEADER_PATH "qlmps")
find_path(PEPS_HEADER_PATH "qlpeps")
find_package(hptt)
find_package(Boost
        1.74.0
        REQUIRED
        COMPONENTS serialization mpi)
find_package(MPI REQUIRED)
#set omp flag
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(OMP_FLAGS -qopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OMP_FLAGS -fopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OMP_FLAGS -fopenmp)
endif ()

#find openmp?
add_compile_options(${MATH_LIB_COMPILE_FLAGS} ${OMP_FLAGS})

include_directories(${MPS_HEADER_PATH}
        ${TENSOR_HEADER_PATH}
        ${PEPS_HEADER_PATH}
        ${hptt_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${MPI_CXX_HEADER_DIR}
)

link_libraries(
        ${hptt_LIBRARY}
        ${MATH_LIB_LINK_FLAGS}
        Boost::serialization
        Boost::mpi
        ${MPI_CXX_LINK_FLAGS}
        ${MPI_mpi_LIBRARY}
)

add_executable(square_simple_update src/square_simple_update.cpp)
add_executable(square_vmc_update src/square_vmc_update.cpp src/myutil.cpp)
add_executable(square_mc_measure src/square_mc_measure.cpp src/myutil.cpp)
add_executable(triangle_simple_update src/triangle_simple_update.cpp)
add_executable(triangle_vmc_update src/triangle_vmc_update.cpp src/myutil.cpp)
add_executable(triangle_mc_measure src/triangle_mc_measure.cpp src/myutil.cpp)
add_executable(triangle_dmrg src/triangle_dmrg.cpp src/myutil.cpp)
add_executable(triangle_dmrg_measure src/triangle_dmrg_measure.cpp src/myutil.cpp)
add_executable(kagome_simple_update src/kagome_simple_update.cpp)
add_executable(kagome_vmc_update src/kagome_vmc_update.cpp src/myutil.cpp)
add_executable(kagome_mc_measure src/kagome_mc_measure.cpp)
add_executable(kagome_dmrg src/kagome_dmrg.cpp src/myutil.cpp)
add_executable(kagome_dmrg_measure src/kagome_dmrg_measure.cpp src/myutil.cpp)
add_executable(transfer_ipeps src/transfer_iPEPS.cpp)
add_executable(increaseD src/tps_increase_D.cpp src/myutil.cpp)
add_executable(extension_tps src/sitps_space_extension.cpp src/myutil.cpp)