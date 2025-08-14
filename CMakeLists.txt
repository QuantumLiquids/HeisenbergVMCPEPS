cmake_minimum_required(VERSION 3.14)
project(HeisenbergVMCPEPS)

set(CMAKE_CXX_STANDARD 20)

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

if (DEFINED U1SYM)
    add_definitions(-DU1SYM=${U1SYM})
    message("USE U1 symmetry code with U1SYM=${U1SYM}.")
endif ()

option(RealCode "tensor element type double or complex" ON)
if (RealCode)
    add_definitions(-DUSE_REAL=1)
else ()
    add_definitions(-DUSE_COMPLEX=1)
endif ()


if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    set(BLA_VENDOR Intel10_64lp_seq)   #FOR VENDOR CHOICE: https://cmake.org/cmake/help/latest/module/FindBLAS.html
    set(BLAS_INCLUDE_DIR "$ENV{MKLROOT}/include")
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64") # (e.g., Apple Silicon)
    add_definitions(-DUSE_OPENBLAS)
    set(BLA_VENDOR OpenBLAS)
    set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas/") #work for homebrew installed openblas in MacOSX
    set(Lapack_ROOT "/opt/homebrew/opt/lapack")
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${OpenBLAS_ROOT} ${Lapack_ROOT}) # For find
    set(BLAS_INCLUDE_DIR "${OpenBLAS_ROOT}/include")
endif ()

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED) # link necessary


find_path(TENSOR_HEADER_PATH "qlten")
find_path(MPS_HEADER_PATH "qlmps")
find_path(PEPS_HEADER_PATH "qlpeps")
find_package(hptt)


include_directories(
        ${MPI_CXX_HEADER_DIR}
        ${BLAS_INCLUDE_DIR}
        ${hptt_INCLUDE_DIR}
        ${TENSOR_HEADER_PATH}
        ${MPS_HEADER_PATH}
        ${PEPS_HEADER_PATH}
)
link_libraries(
        ${hptt_LIBRARY}
        BLAS::BLAS
        LAPACK::LAPACK
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX
)

add_library(myutil_obj OBJECT src/myutil.cpp)

add_executable(square_simple_update src/square_simple_update.cpp)
add_executable(square_vmc_update src/square_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(square_vmc_update_migrated src/square_vmc_update_migrated.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(square_mc_measure src/square_mc_measure.cpp $<TARGET_OBJECTS:myutil_obj>)

add_executable(square_planar_simple_update src/square_planar_simple_update.cpp)
add_executable(square_planar_vmc_update src/square_planar_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)


add_executable(triangle_simple_update src/triangle_simple_update.cpp)
add_executable(triangle_loop_update src/triangle_loop_update.cpp)
add_executable(triangle_vmc_update src/triangle_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_vmc_update_no_u1 src/triangle_vmc_update_no_u1.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_mc_measure src/triangle_mc_measure.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_vmps src_dmrg/triangle_vmps.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_dmrg src_dmrg/triangle_dmrg.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_dmrg_measure src_dmrg/triangle_dmrg_measure.cpp $<TARGET_OBJECTS:myutil_obj>)
add_executable(triangle_dmrg_measure_sf src_dmrg/triangle_dmrg_measure_sf.cpp $<TARGET_OBJECTS:myutil_obj>)


#add_executable(kagome_simple_update src/kagome_simple_update.cpp)
#add_executable(kagome_vmc_update src/kagome_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(kagome_mc_measure src/kagome_mc_measure.cpp)
#add_executable(kagome_dmrg src_dmrg/kagome_dmrg.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(kagome_dmrg_measure src_dmrg/kagome_dmrg_measure.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(transfer_ipeps src/transfer_iPEPS.cpp)

# Utility tools
if (NOT DEFINED U1SYM)
    add_executable(increaseD src/tps_increase_D.cpp $<TARGET_OBJECTS:myutil_obj>)
    add_executable(extend src/sitps_space_extension.cpp $<TARGET_OBJECTS:myutil_obj>)
endif ()
add_executable(zero_evolve src/square_peps_zero_evolve.cpp)
#add_executable(fluctuation src_vmc_boost/fluctuation_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)