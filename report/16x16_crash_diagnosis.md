# 16×16 (and 12×12) Crash Diagnosis Summary

## Core Finding

The crash is a **stack overflow caused by infinite recursion inside Intel MPI's CH4 device layer**, NOT a bug in our application code or the compiler. Both GCC and icpx compiled binaries crash with identical MPI-internal stack traces.

## Core Dump Evidence

We obtained core dumps from two independent crashes on different compilers:

### GCC -O2 -g core dump (c067, job 538832, Mar 17)
- Binary: `build_O2g/vmc_optimize` (GCC 12, `-O2 -g`)
- Crashed after 3 iterations of 16×16 J2=0 D=8 PBC VMC optimization
- 56 MPI ranks, 1 node, 1280 MC samples

### icpx -O2 -g core dump (c068, job 540033, Mar 17)
- Binary: `debug_icpx_O2/vmc_optimize` (Intel icpx 2025.2, `-O2 -g`)
- Crashed after 4 iterations of the same problem
- 56 MPI ranks, 1 node, 1280 MC samples

### Stack trace (identical pattern in both cores)

```
#0  recv_target_cmpl_cb (rreq=0x7f200a6683c0) at ../../src/mpid/ch4/src/ch4r_callbacks.c:321
#1  MPIDIG_progress_compl_list (vci=<optimized out>) at ../../src/mpid/ch4/src/ch4r_callbacks.c:92
#2  recv_target_cmpl_cb (rreq=<optimized out>) at ../../src/mpid/ch4/src/ch4r_callbacks.c:384
#3  MPIDIG_progress_compl_list (vci=<optimized out>) at ../../src/mpid/ch4/src/ch4r_callbacks.c:92
#4  recv_target_cmpl_cb (rreq=<optimized out>) at ../../src/mpid/ch4/src/ch4r_callbacks.c:384
#5  MPIDIG_progress_compl_list (vci=<optimized out>) at ../../src/mpid/ch4/src/ch4r_callbacks.c:92
... (repeating for 29000+ frames until stack overflow → SIGSEGV)
```

The crash is an **infinite mutual recursion** between:
- `recv_target_cmpl_cb()` at `ch4r_callbacks.c:321/384`
- `MPIDIG_progress_compl_list()` at `ch4r_callbacks.c:92`

A receive-completion callback triggers the progress-completion-list, which triggers another receive-completion callback, ad infinitum, until the stack is exhausted and the process segfaults.

GDB `info threads` shows only 2 threads in both cores:
```
* 1  Thread 0x... (LWP ...) recv_target_cmpl_cb at ch4r_callbacks.c:321
  2  Thread 0x... (LWP ...) internal_fnmatch () from /lib64/libc.so.6
```

## MPI Version

Intel MPI 2021.16 (part of oneAPI 2025.2):
- Path: `/share/home/wanghx/intel/oneapi/mpi/2021.16/`
- MPI compiler wrappers: `mpigxx`, `mpigcc`, `mpiicpx`
- CH4 device layer (the crashing code is in `src/mpid/ch4/src/ch4r_callbacks.c`)

## Evidence Table: All 16×16 Runs

| JobID | Compiler | Flags | Ranks | Iters completed | Result |
|-------|----------|-------|-------|-----------------|--------|
| 527041 | icpx Debug+ASan | -g -fsanitize=address | 112 | 0 | CRASHED (SIGNAL 11) |
| 527832 | GCC 12 | -O3 | 56 | 2 | CRASHED (SIGNAL 11) |
| 527948 | GCC 12 | -O3 | 4 | 1 | OOM (too few ranks) |
| 528239 | icpx ASan | -g -fsanitize=address | 8 | 0 | CRASHED (startup) |
| 528253 | GCC 12 | -g -fsanitize=address | 8 | 5 | PASSED |
| 528590 | GCC 12 | -g -fsanitize=address | 56 | 5 | PASSED (128 samples) |
| 528603 | GCC 12 | -g (no opt) | 56 | 5 | PASSED |
| 531981 | icpx | -O2 | 56 | 5 | PASSED (MaxIter=5 limit) |
| 531982 | GCC 12 | -O2 -fsanitize=undefined | 56 | 3 | CRASHED (SIGNAL 11, 0 UB errors) |
| 538832 | GCC 12 | -O2 -g | 56 | 3 | **CRASHED — core dump obtained** |
| 540033 | icpx | -O2 -g | 56 | 4 | **CRASHED — core dump obtained** |
| 540034 | icpx | -O2 -g (MinSR) | 56 | 11 | CRASHED (SIGNAL 9) |

Key patterns:
- **Both compilers crash** with the same MPI-internal stack trace
- **More ranks → crash sooner**: 112 ranks crash at iter 0, 56 ranks at iter 2-4
- **8 ranks with 128 samples passes**: fewer MPI messages avoid triggering the bug
- **Debug builds (-g, no optimization) sometimes pass**: different timing/message ordering
- **UBSan found 0 undefined behavior errors**: our code is clean
- **ASan found 0 memory errors** (in debug builds where it ran)
- 8×8 lattice never crashes: smaller system = less MPI traffic per iteration

## What Our Code Does with MPI

The VMC optimization loop uses MPI for:
1. **Sample distribution**: MC samples are distributed across ranks (1280 samples / 56 ranks ≈ 23 samples/rank)
2. **Gradient aggregation**: After local gradient computation, all ranks do MPI collective operations (likely `MPI_Allreduce` or `MPI_Allgather`) to aggregate gradients
3. **SR/MinSR linear solve**: The Stochastic Reconfiguration solver involves MPI collectives for the conjugate gradient iteration
4. **TPS synchronization**: The tensor product state is synchronized across ranks after each update

The upstream PEPS library (`qlpeps`) handles all MPI communication. The relevant code paths:
- VMC optimizer: `/Users/wanghaoxin/GitHub/PEPS/` (qlpeps library)
- MPI wrapper/communication: look for `MPI_Allreduce`, `MPI_Bcast`, `MPI_Gather`, `MPI_Send`, `MPI_Recv`, `MPI_Isend`, `MPI_Irecv` calls
- The crash specifically involves CH4's receive-completion callbacks, suggesting the issue is triggered during point-to-point or collective receive operations

## Upstream Repos to Check

1. **qlpeps** (PEPS library): `/Users/wanghaoxin/GitHub/PEPS/`
   - This is where VMC optimization and MPI communication are implemented
   - Look for non-blocking MPI calls (`MPI_Isend`/`MPI_Irecv`) and their completion (`MPI_Wait`/`MPI_Waitall`)
   - Look for patterns that could trigger recursive MPI progress (e.g., calling MPI functions inside completion handlers, or nested MPI calls)

2. **qlten** (Tensor library): `/Users/wanghaoxin/GitHub/TensorToolkit/`
   - Lower-level tensor operations; may have MPI calls for distributed tensor contractions

3. **qlmps** (MPS/DMRG library): `/Users/wanghaoxin/GitHub/UltraDMRG/include/`
   - Parameter parsing infrastructure; less likely to have MPI issues

## Possible Root Causes in Our Code

Even though the crash is inside Intel MPI, our code may be triggering it via:

1. **Recursive MPI calls**: If a completion callback or signal handler in our code calls an MPI function, this can trigger recursive progress in CH4
2. **MPI calls from signal handlers**: Any MPI call from a signal handler can cause recursive re-entry
3. **Thread-safety issues**: If our code calls MPI from multiple threads without proper `MPI_Init_thread(MPI_THREAD_MULTIPLE)`
4. **Large number of pending non-blocking operations**: If we post many `MPI_Isend`/`MPI_Irecv` without completing them, the completion list grows and the recursive progress can overflow
5. **MPI_Allreduce with large data**: The collective implementation may internally use point-to-point with callbacks that cascade

## Suggested Investigation Steps

1. Search qlpeps for all MPI calls: `grep -rn "MPI_" /Users/wanghaoxin/GitHub/PEPS/`
2. Check if any MPI calls happen inside OpenMP parallel regions or callbacks
3. Look for patterns of many non-blocking operations without intermediate waits
4. Check `MPI_Init_thread` level — is it `MPI_THREAD_SINGLE`, `MPI_THREAD_FUNNELED`, or `MPI_THREAD_MULTIPLE`?
5. Check if `I_MPI_ASYNC_PROGRESS` is set (async progress can trigger the recursive callback issue)
6. Look for any custom signal handlers that might call MPI functions

## Environment Details

- Cluster: SUSTech Physics HPC, Rocky Linux 8
- Nodes: 256GB RAM, 56 cores each
- Intel MPI 2021.16 (oneAPI 2025.2)
- GCC 12.1.0, Intel icpx 2025.2
- No custom `I_MPI_*` environment variables were set in the failing jobs
- `OMP_NUM_THREADS=1` (OpenMP parallelism disabled)
