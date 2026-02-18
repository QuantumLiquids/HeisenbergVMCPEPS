---
name: vmcpeps-cluster-workflow
description: "Use when the user asks to run calculations, submit jobs, check
job status, fetch results, or manage cases on the susphy HPC cluster. Covers
the full pipeline: simple update, tile, VMC optimize, MC measure."
---

# VMC-PEPS Cluster Workflow

Use this skill when setting up, submitting, monitoring, or fetching results
for VMC-PEPS calculations on the susphy cluster.

## Pipeline Overview

Three slurm jobs with dependency chaining:

1. **Stage 1 (SU + tile):** Single-node. Runs `simple_update` on unit cell
   with tau-schedule + advanced-stop, then `sitps_tile` to expand to full
   lattice.
2. **Stage 2 (VMC):** MPI job. Runs `vmc_optimize` on full lattice.
   `--dependency=afterok:stage1`.
3. **Stage 3 (Measure):** MPI job. Runs `mc_measure` on optimized state.
   `--dependency=afterok:stage2`.

## CLI Reference

All commands via `python3 scripts/cluster_run.py`:

### Generate a new case

```bash
python3 scripts/cluster_run.py new-case \
  --lattice 8x8 --D 8 --J2 0.0 --bc pbc \
  --unit-cell 2x2 \
  --partition 256G56c --ntasks 56
```

Key optional flags: `--optimizer SR`, `--lr 0.1`, `--iterations 30`,
`--taus "0.5,0.2,0.1,0.05,0.02"`, `--override key=value`.

### Submit to cluster

```bash
python3 scripts/cluster_run.py submit \
  --case run/8x8J2=0D8/<timestamp>
```

### Check status

```bash
python3 scripts/cluster_run.py status
```

Or directly: `ssh susphy "squeue -u wanghx"`

### Fetch results

```bash
python3 scripts/cluster_run.py fetch \
  --case run/8x8J2=0D8/<timestamp> --output-dir data/
```

## Pre-Submit Checklist

Before submitting any job:

1. Show the full parameter summary to the user
2. Confirm the user approves the parameters (especially D, lattice size,
   optimizer, iteration count)
3. Verify SSH connectivity: `ssh susphy "hostname"`
4. Check no conflicting jobs: `ssh susphy "squeue -u wanghx"`

## Parameter Defaults

- MC samples (VMC): `max(10000, min(20000, Lx*Ly*200))`, rounded to ntasks
- MC samples (measure): 2x VMC samples
- TRGDmin = D, TRGDmax = 2*D (PBC only)
- Dbmps_min = D, Dbmps_max = 2*D (OBC only)
- Tau schedule: `0.5, 0.2, 0.1, 0.05, 0.02`
- Advanced stop: energy_tol=1e-6, lambda_tol=1e-5, patience=3

## Cluster Operational Notes (susphy)

### Partition selection
- Use `256G56c` for production runs. The 24c partitions have ABI
  compatibility issues with login-node builds.

### Environment loading
Slurm scripts must source the environment with nounset temporarily disabled:

```bash
set +u; set +x
source /share/home/wanghx/myenv.sh
set -x; set -u
```

### Failure signatures
- `FAILED 0:00 with empty output on 64G24c`: partition unusable, switch to
  256G56c.
- `libmkl_* not found` / `GLIBCXX_* not found`: env not loaded in job script.
- `BASHRCSOURCED unbound`: don't source .bashrc under set -u.

### Expected success checkpoints
- `Simple Update completed.`
- `SplitIndexTPS saved to: tpsfinal`
- `Optimization completed!`
- `Rank 0: statistic data finished.`

## Result Files

After a successful run, key output files are:

- `vmc/energy/energy_trajectory.csv` -- columns: iteration, energy,
  energy_error, gradient_norm
- `measure/stats/energy.csv` -- final energy: index, mean, stderr
- `measure/stats/SzSz_all2all.csv` -- spin-spin correlations
- `measure/stats/bond_energy_{h,v}_{mean,stderr}.csv` -- bond energies

## Directory Structure

```
run/<case>/<timestamp>/
  params/         -- frozen input JSONs
  stage1_su_tile.slurm
  stage2_vmc.slurm
  stage3_measure.slurm
  su_unit/        -- simple update working dir
  tpsfinal_tiled/ -- tiled TPS
  vmc/            -- VMC working dir (energy/, tpsfinal/, checkpoints/)
  measure/        -- measurement working dir (stats/, tpsfinal/)
  logs/           -- slurm stdout/stderr
  meta.json       -- run metadata + job IDs
```
