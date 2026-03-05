---
name: research-notes
description: "Research log for HeisenbergVMCPEPS calculations. Use when discussing
  results, comparing energies, planning new runs, or summarizing progress. Contains
  accumulated findings, energy tables, and physics insights."
---

# Research Notes: Square Heisenberg PEPS+VMC

## Current Energy Summary

### J2 = 0 (Pure Heisenberg)

| System | D | Source | Method | E/site | Error | Status | Notes |
|--------|---|--------|--------|--------|-------|--------|-------|
| 8x8 | 5 | Yubin cluster update | Verify (measure only) | -0.67181 | 0.00003 | Done | Baseline |
| 8x8 | 5 | Yubin cluster update | VMC 30 iters | -0.67317 | 0.000033 | Done | +0.2% improvement |
| 8x8 | 8 | Yubin cluster update | Verify (measure only) | -0.67187 | -- | Done | Similar to D=5 Yubin |
| 8x8 | 8 | Yubin cluster update | VMC (running) | -- | -- | Running | Job 523623 |
| 12x12 | 8 | Yubin cluster update | VMC (running) | -- | -- | Running | Job 523634 |
| 16x16 | 8 | Yubin cluster update | VMC 2-node | -- | -- | Running | Job 524117. Single-node OOM'd (220 GB) |
| 8x8 | -- | -- | QMC benchmark | -0.673487 | -- | Ref | |
| 12x12 | -- | -- | QMC benchmark | -0.670685 | -- | Ref | |
| 16x16 | -- | -- | QMC benchmark | -0.669976 | -- | Ref | |

### J2 = 0.5 (Frustrated)

| System | D | Source | Method | E/site | Error | Status | Notes |
|--------|---|--------|--------|--------|-------|--------|-------|
| 8x8 | 8 | Our SU | SU initial (VMC iter 0) | -0.496 | -- | Done | From 20260220 run |
| 8x8 | 8 | Yubin cluster update | Verify (measure only) | -0.49648 | -- | Done | Only marginally better than SU |
| 8x8 | 5 | Yubin cluster update | Verify (measure only) | -0.49490 | 0.00006 | Done | Baseline |
| 8x8 | 5 | Yubin cluster update | VMC 30 iters | -0.49767 | -- | Done | +0.6% improvement |
| 8x8 | 8 | Yubin cluster update | VMC (running) | -- | -- | Running | Job 524096 |
| 8x8 | 8 | Our SU (iter 4) | VMC continue (running) | -- | -- | Running | Job 524679. Resumed from iter 4 |
| 16x16 | 8 | Yubin cluster update | VMC 2-node | -- | -- | Running | Job 524119 |

## Key Findings

### SU vs Yubin tensor quality

- **J2=0.5 D=8**: Yubin's cluster-update tensor (E/site = -0.49648) is only marginally
  better than our simple update initial state (E/site ~ -0.496). The gap is
  ~0.001 in E/site. This suggests cluster-update optimization at J2=0.5 converges to
  a state close to what SU already achieves.
- **J2=0 D=5**: Yubin's tensor (-0.67181) is below our SU starting points, and
  VMC further improves it to -0.67317 (within 0.05% of QMC at -0.673487).
- **J2=0 D=8 vs D=5**: D=8 Yubin verify (-0.67187) is very close to D=5
  Yubin verify (-0.67181). The D increase alone doesn't help much without VMC.

### VMC improvement over baselines

- J2=0 D=5: VMC improves Yubin baseline by ~0.2% (from -0.67181 to -0.67317)
- J2=0.5 D=5: VMC improves Yubin baseline by ~0.6% (from -0.49490 to -0.49767)
- Frustrated (J2=0.5) models benefit more from VMC optimization

### Memory scaling for SR optimizer

The dominant memory cost is O* (log-derivative) sample storage:

```
Memory per rank ~ N_samples_per_rank * N_sites * d * D^4 * 8 bytes
```

Actual measurements from sacct/sstat (56 ranks per node):

| System | D | Stage | MaxRSS (node) | Per rank | Source |
|--------|---|-------|---------------|----------|--------|
| 8x8 | 5 | VMC | 20.7 GB | 370 MB | sacct (job 522522) |
| 8x8 | 8 | Measure | 13.6 GB | 243 MB | sacct (job 524095) |
| 8x8 | 8 | VMC | 70.7 GB | 1.26 GB | sstat (job 523623) |
| 16x16 | 8 | VMC (OOM) | 220.5 GB | 3.9 GB | sacct (job 523636) |

Memory scaling observations:
- D=5→D=8 VMC: 20.7→70.7 GB (3.4x). Pure D^4 ratio would be 6.6x.
  Memory doesn't scale purely as D^4 — there are fixed overheads.
- 8x8→16x16 D=8 VMC: 70.7→220+ GB (3.1x). Sites ratio is 4x.
  Sublinear in sites, likely because O* samples dominate only above
  a threshold; TRG contraction adds a smaller per-site cost.
- Measure uses much less memory than VMC (~13.6 vs 70.7 GB for D=8 8x8)
  because measure doesn't store O* samples.

**Conclusion**: 16x16 D=8 requires 2 nodes (512 GB). For D>8 at 16x16, need
further upstream optimization (streaming S*v, minibatch SR).

### Loop update

- Attempted for J2=0 D=5 on 8x8. Energy went up compared to simple update.
  Energy oscillation from normalization factor is acceptable, but the overall
  result was worse. Root cause unclear. Not pursued further for now.

## Benchmark References

QMC benchmarks stored in `data/benchmarks/qmc_square_heisenberg.json`.

Currently available:
- **J2=0**: L=6 (-0.6788), L=8 (-0.673487), L=12 (-0.670685), L=16 (-0.669976)
- **J2=0.5**: No QMC benchmarks yet (frustrated model, sign problem)

## Run Directory Index

All runs at `run/{physics_case}/{timestamp_suffix}/`. Each has `run_history.json`.

### Completed runs worth referencing

| Run directory | Purpose | Key result |
|---------------|---------|------------|
| 8x8J2=0D5/20260220_yubin_verify | D=5 J2=0 Yubin baseline | E/site = -0.67181 |
| 8x8J2=0D5/20260220_201619_yubin_vmc | D=5 J2=0 VMC | E/site = -0.67317 |
| 8x8J2=0.5D5/20260220_yubin_verify | D=5 J2=0.5 Yubin baseline | E/site = -0.49490 |
| 8x8J2=0.5D5/20260223_yubin_vmc | D=5 J2=0.5 VMC | E/site = -0.49767 |
| 8x8J2=0D8/20260227_yubin_verify | D=8 J2=0 Yubin baseline | E/site = -0.67187 |
| 8x8J2=0.5D8/20260301_yubin_verify | D=8 J2=0.5 Yubin baseline | E/site = -0.49648 |
| 8x8J2=0.5D8/20260220_013704_su_vmc_measure | D=8 J2=0.5 SU state | E/site ~ -0.496 (iter 0) |
| 8x8J2=0.5D8/20260302_su_vmc_continue | D=8 J2=0.5 SU VMC continue | Running (job 524679) |

### Failed runs with lessons

| Run directory | What happened | Lesson |
|---------------|--------------|--------|
| 8x8J2=0D5/20260220_011019_su_loopupdate | Loop update raised energy above SU | Root cause unclear |
| 16x16J2=0D8/20260227_yubin_vmc | OOM at 220 GB | 16x16 D=8 needs 2 nodes |
| 8x8J2=0D8/20260217_211901_su_vmc | Timeout, 28 ranks | Always use 56 ranks per node |

## Open Questions

- Will D=8 VMC significantly improve over D=5 VMC for J2=0? (D=8 baseline barely
  differs from D=5 baseline, but VMC may find better states with more parameters)
- What is the finite-size scaling of PEPS+VMC energy toward thermodynamic limit?
