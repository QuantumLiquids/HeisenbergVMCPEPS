---
name: research-notes
description: "Research log for HeisenbergVMCPEPS calculations. Use when discussing
  results, comparing energies, planning new runs, or summarizing progress. Contains
  accumulated findings, energy tables, and physics insights."
---

# Research Notes: Square Heisenberg PEPS+VMC

## Production Run Strategy (established 2026-03-19)

- **Initial TPS**: Always Yubin's 2×2 cluster-update tensor from `yubin_tensor/` tiled via `sitps_tile`
- **WarmUp**: Always set `WarmUp=100` (program skips if configs exist)
- **Checkpoint**: Every 1 iteration (`CheckpointEveryNSteps=1`)
- **Resume**: On crash, resubmit from last tpsfinal + configs
- **Binary**: `build_icpx_prod/vmc_optimize` (icpx 2025.3.2, Intel MPI 2021.17, `-O2 -g`)
- **Standard params**: 1280 samples, lr=0.1, SRDiagShift=0.01, TRGDmin=8, TRGDmax=16

## Current Production Runs (2026-03-19)

### J2 = 0

| JobID | System | Optimizer | Initial TPS | Configs | Status |
|-------|--------|-----------|-------------|---------|--------|
| 554682 | 16×16 | MinSR cont | from MinSR iter 10 | yes | RUNNING |
| 554726 | 16×16 | SR fresh | Yubin CU tiled | yes (from MinSR) | RUNNING |
| 554840 | 12×12 | SR fresh | Yubin CU tiled | yes (from 540093) | RUNNING |
| 554841 | 12×12 | MinSR fresh | Yubin CU tiled | yes (from 540093) | RUNNING |

### J2 = 0.5

| JobID | System | Optimizer | Initial TPS | Configs | Status |
|-------|--------|-----------|-------------|---------|--------|
| 554800 | 12×12 | SR | Yubin CU tiled | no (WarmUp=100) | RUNNING |
| 554801 | 12×12 | MinSR | Yubin CU tiled | no (WarmUp=100) | RUNNING |
| 554802 | 16×16 | SR | Yubin CU tiled | no (WarmUp=100) | RUNNING |
| 554803 | 16×16 | MinSR | Yubin CU tiled | no (WarmUp=100) | RUNNING |

## Completed 8×8 Results

### J2 = 0

| System | D | Source | Method | E/site | Error | Run Dir | Job | Notes |
|--------|---|--------|--------|--------|-------|---------|-----|-------|
| 8x8 | 5 | Yubin CU | Verify (measure) | -0.67181 | 0.00003 | `8x8J2=0D5/20260220_yubin_verify` | - | Baseline |
| 8x8 | 5 | Yubin CU | VMC 30 iters | -0.67317 | 0.000033 | `8x8J2=0D5/20260220_201619_yubin_vmc` | - | +0.2% vs verify |
| 8x8 | 8 | Yubin CU | Verify (measure) | -0.67187 | -- | `8x8J2=0D8/20260227_yubin_verify` | - | Similar to D=5 |
| 8x8 | 8 | Yubin CU | SR 30 iter | -0.67342 | -- | `8x8J2=0D8/20260303_yubin_vmc_stepselector` | 525891 | |
| 8x8 | 8 | Yubin CU | StepSelector 60 iter | **-0.67349** | -- | `8x8J2=0D8/20260303_yubin_vmc_stepselector` | 525891 | 3-phase lr: 0.3→0.15→0.075. Gap to QMC < 10⁻⁵ |
| 8x8 | 8 | Yubin CU | MinSR lr=0.1 60 iter | -0.67340 | -- | `8x8J2=0D8/` (MinSR run) | - | Validates MinSR. lr=0.3 diverges |
| 8x8 | -- | -- | QMC benchmark | -0.673487 | -- | -- | -- | Reference |
| 8x8 | 5 | Yubin | GO (Yubin's method) | -0.6733 | -- | -- | -- | ~28 steps |

### J2 = 0.5

| System | D | Source | Method | E/site | Error | Run Dir | Job | Notes |
|--------|---|--------|--------|--------|-------|---------|-----|-------|
| 8x8 | 8 | Our SU | SU initial (VMC iter 0) | -0.496 | -- | `8x8J2=0.5D8/20260220_013704_su_vmc_measure` | - | From SU pipeline |
| 8x8 | 8 | Yubin CU | Verify (measure) | -0.49648 | -- | `8x8J2=0.5D8/20260301_yubin_verify` | 524095 | Marginally better than SU |
| 8x8 | 5 | Yubin CU | Verify (measure) | -0.49490 | 0.00006 | `8x8J2=0.5D5/20260220_yubin_verify` | - | Baseline |
| 8x8 | 5 | Yubin CU | VMC 30 iters | -0.49767 | -- | `8x8J2=0.5D5/20260223_yubin_vmc` | - | +0.6% vs verify |
| 8x8 | 8 | Yubin CU | SR 30 iter | -0.49814 | -- | `8x8J2=0.5D8/20260302_su_vmc_continue` | 524679 | |
| 8x8 | 8 | Yubin CU | StepSelector 30 iter | **-0.49840** | -- | `8x8J2=0.5D8/20260307_yubin_vmc_stepselector` | 527811 | lr stayed at 0.3 (AutoSelector) |
| 8x8 | 8 | Yubin | GO (Yubin's method) | -0.4978 | -- | -- | -- | ~15 steps |

### Key 8×8 Findings
- **SU vs Yubin CU (J2=0.5 D=8)**: Yubin CU (-0.49648) only marginally better than SU (-0.496). Gap ~0.001.
- **J2=0 D=8 vs D=5**: D=8 verify (-0.67187) ≈ D=5 verify (-0.67181). D increase alone doesn't help without VMC.
- **VMC D=5 J2=0**: 30 iters improves from -0.67181 to -0.67317 (+0.2%), within 0.05% of QMC.
- **StepSelector**: Reaches QMC-level on J2=0 (gap < 10⁻⁵). Beats Yubin on J2=0.5 (-0.4984 vs -0.4978).

## Completed 16×16 Results (new MPI, 2026-03-24)

### J2 = 0
| System | D | Source | Method | E/site | Run Dir | Job | Notes |
|--------|---|--------|--------|--------|---------|-----|-------|
| 16x16 | 8 | Yubin CU | SR 60 iter | **-0.66967** | `16x16J2=0D8/20260319_icpx_sr_fresh` | 554726 | 4d 11h. Best at iter 45. Gap to QMC: 4.6×10⁻⁴ |
| 16x16 | 8 | Yubin CU | MinSR 60 iter | **-0.66967** | `16x16J2=0D8/20260319_icpx_minsr_fresh` | 554923 | 4d 14h. Best at iter 34. Matches SR exactly |
| 16x16 | -- | -- | QMC | -0.669976 | -- | -- | Reference |
| 16x16 | 5 | Yubin | GO | -0.6696 | -- | -- | Yubin's result |

Both SR and MinSR surpass Yubin D=5 (-0.6696) by 0.01%. First successful 60-iter 16×16 runs after MPI bug fix.

## In-Progress 12×12 and 16×16 J2=0.5

| System | J2 | Method | Iters | Best E/site | Run Dir | Job |
|--------|-----|--------|-------|-------------|---------|-----|
| 12x12 | 0 | SR | 12/60 | -0.67007 | `12x12J2=0D8/20260319_icpx_sr_fresh` | 554840 |
| 12x12 | 0 | MinSR | 12/60 | -0.67009 | `12x12J2=0D8/20260319_icpx_minsr_fresh` | 554841 |
| 12x12 | 0.5 | SR | 7/60 | -0.49585 | `12x12J2=0.5D8/20260319_icpx_sr_60iter` | 554800 |
| 12x12 | 0.5 | MinSR | 7/60 | -0.49585 | `12x12J2=0.5D8/20260319_icpx_minsr_60iter` | 554801 |
| 16x16 | 0.5 | SR | 51/60 | -0.49588 | `16x16J2=0.5D8/20260319_icpx_sr_60iter` | 554802 |
| 16x16 | 0.5 | MinSR | 52/60 | -0.49573 | `16x16J2=0.5D8/20260319_icpx_minsr_60iter` | 554803 |

## Old Partial 12×12 / 16×16 Results (pre-MPI-fix, for reference)

### 12×12 J2=0 (crashed after 7 iters, GCC)
- iter 0: -0.67106 (biased low, WarmUp=10), iter 6: -0.67018
- QMC: -0.670685, Yubin D=5: -0.6704

### 16×16 J2=0 (multiple runs, all crashed due to MPI bug)
- Best SR (531981+540033, 9 combined iters): -0.66929
- Best MinSR (540034, 11 iters): -0.66954
- QMC: -0.669976, Yubin D=5: -0.6696

## 16×16 Crash Diagnosis

**Root cause**: Infinite recursion in Intel MPI CH4 layer (`recv_target_cmpl_cb` ↔ `MPIDIG_progress_compl_list`), causing stack overflow → SIGSEGV. NOT compiler-specific, NOT our code. See `report/16x16_crash_diagnosis.md`.

**Fix**: Upstream qlpeps MPI usage improved (two-phase broadcast, packed wire formats) + Intel MPI updated from 2021.16 to 2021.17. New production runs use updated code.

## Key Findings

- **SU vs Yubin CU tensor**: Yubin's cluster-update tensor is better initial state than simple update. Always use it.
- **WarmUp**: Must be ≥100 for fresh starts. WarmUp=0 or WarmUp=10 gives biased iter 0 energy.
- **MinSR**: Works with lr=0.1 (lr=0.3 diverges). Solves N_samples² system vs SR's N_params². Critical for large D.
- **StepSelector**: Effective on 8×8 — reaches QMC-level accuracy. Not yet tested on larger systems.
- **16×16 crash**: Was MPI bug, not code or compiler. Fixed by upstream MPI improvements.

## Reference Values

| System | QMC (J2=0) | Yubin D=5 (J2=0) | Yubin D=8 (J2=0.5) |
|--------|------------|-------------------|---------------------|
| 8×8 | -0.673487 | -0.6733 | -0.4978 |
| 12×12 | -0.670685 | -0.6704 | ~-0.497 |
| 16×16 | -0.669976 | -0.6696 | ~-0.496 |
