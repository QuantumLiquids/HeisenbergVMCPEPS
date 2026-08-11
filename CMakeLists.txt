cmake_minimum_required(VERSION 3.27)
project(HeisenbergVMCPEPS LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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


find_package(PEPS 0.1 CONFIG REQUIRED)

# PEPS::PEPS carries the complete public dependency chain, including
# UltraDMRG, TensorToolkit, MPI, the numeric backend, HPTT, and optional
# ScaLAPACK support selected when the PEPS package was built.
function(add_heisenberg_executable target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE PEPS::PEPS)
endfunction()

add_library(myutil_obj OBJECT src/myutil.cpp)
target_link_libraries(myutil_obj PRIVATE PEPS::PEPS)

add_heisenberg_executable(simple_update src/simple_update.cpp)
add_heisenberg_executable(loop_update src/loop_update.cpp)
add_heisenberg_executable(vmc_optimize src/vmc_optimize.cpp)
add_heisenberg_executable(mc_measure src/mc_measure.cpp)

## Removed deprecated per-model entries; use unified simple_update/vmc_optimize/mc_measure


## NOTE: Deprecated; triangle now handled by unified simple_update (ModelType="TriangleHeisenberg").
option(BUILD_TRIANGLE_LOOP_UPDATE "Build triangle loop update tool (legacy)" OFF)
if (BUILD_TRIANGLE_LOOP_UPDATE)
    add_heisenberg_executable(triangle_loop_update src/triangle_loop_update.cpp)
endif ()

option(BUILD_TRIANGLE_VMC_NO_U1 "Build triangle VMC no-U1 tool (legacy API)" OFF)
if (BUILD_TRIANGLE_VMC_NO_U1)
    add_heisenberg_executable(triangle_vmc_update_no_u1 src/triangle_vmc_update_no_u1.cpp $<TARGET_OBJECTS:myutil_obj>)
endif ()
add_heisenberg_executable(triangle_vmps src_dmrg/triangle_vmps.cpp $<TARGET_OBJECTS:myutil_obj>)
add_heisenberg_executable(triangle_dmrg src_dmrg/triangle_dmrg.cpp $<TARGET_OBJECTS:myutil_obj>)
add_heisenberg_executable(triangle_dmrg_measure src_dmrg/triangle_dmrg_measure.cpp $<TARGET_OBJECTS:myutil_obj>)
add_heisenberg_executable(triangle_dmrg_measure_sf src_dmrg/triangle_dmrg_measure_sf.cpp $<TARGET_OBJECTS:myutil_obj>)


#add_executable(kagome_simple_update src/kagome_simple_update.cpp)
#add_executable(kagome_vmc_update src/kagome_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(kagome_mc_measure src/kagome_mc_measure.cpp)
#add_executable(kagome_dmrg src_dmrg/kagome_dmrg.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(kagome_dmrg_measure src_dmrg/kagome_dmrg_measure.cpp $<TARGET_OBJECTS:myutil_obj>)
#add_executable(transfer_ipeps src/transfer_iPEPS.cpp)

# Utility tools
if (NOT DEFINED U1SYM)
    add_heisenberg_executable(increaseD src/tps_increase_D.cpp $<TARGET_OBJECTS:myutil_obj>)
endif ()
add_heisenberg_executable(sitps_tile src/sitps_tile.cpp $<TARGET_OBJECTS:myutil_obj>)
option(BUILD_ZERO_EVOLVE "Build zero-evolve tool (legacy)" OFF)
if (BUILD_ZERO_EVOLVE)
    add_heisenberg_executable(zero_evolve src/square_peps_zero_evolve.cpp)
endif ()
#add_executable(fluctuation src_vmc_boost/fluctuation_vmc_update.cpp $<TARGET_OBJECTS:myutil_obj>)

# Optional (temporary) tools — build only when needed
option(BUILD_CONVERTER "Build TPS->SplitIndexTPS converter (temporary tool)" OFF)
if (BUILD_CONVERTER)
    add_heisenberg_executable(convert_tps_to_sitps src/convert_tps_to_sitps.cpp $<TARGET_OBJECTS:myutil_obj>)
endif ()
option(BUILD_TXT_TO_TPS_SITPS "Build TXT->TPS->SplitIndexTPS converter (temporary tool)" OFF)
if (BUILD_TXT_TO_TPS_SITPS)
    add_heisenberg_executable(txt_to_tps_sitps src/txt_to_tps_sitps.cpp)
endif ()
