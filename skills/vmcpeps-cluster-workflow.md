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
5. Create or update `run_history.json` with job entry (`status: "RUNNING"`)

## Parameter Defaults

- MC samples (VMC): `max(10000, min(20000, Lx*Ly*200))`, rounded to ntasks
- MC samples (measure): 2x VMC samples
- TRGDmin = D, TRGDmax = 2*D (PBC only)
- Dbmps_min = D, Dbmps_max = 2*D (OBC only)
- Tau schedule: `0.5, 0.2, 0.1, 0.05, 0.02`
- Advanced stop: energy_tol=1e-6, lambda_tol=1e-5, patience=3

## Cluster Operational Notes (susphy)

### Core configuration (256G56c)
- **56 cores per node**, 256 GB RAM.
- VMC and measure (MPI stages): set `--ntasks` to `56 * N` where N is
  the number of nodes (`--nodes N`), with `OMP_NUM_THREADS=1`. Each MPI
  rank gets one core. Using fewer ranks wastes cores.
- Simple update (stage 1) is single-process; core count irrelevant.

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
run/<physics_case>/<timestamp_suffix>/
  run_history.json -- structured job history (see Run History Protocol)
  params/          -- frozen input JSONs
  stage1_su_tile.slurm
  stage2_vmc.slurm
  stage3_measure.slurm
  su_unit/         -- simple update working dir
  tpsfinal_tiled/  -- tiled TPS
  vmc/             -- VMC working dir (energy/, tpsfinal/, checkpoints/)
  measure/         -- measurement working dir (stats/, tpsfinal/)
  logs/            -- slurm stdout/stderr
```

Top-level directories are named by physics parameters (e.g., `8x8J2=0D5/`).
Date subdirectories have descriptive suffixes (e.g., `20260220_yubin_verify/`,
`20260220_201619_yubin_vmc/`).

## Run Journal Protocol

A persistent `run_journal.jsonl` (gitignored) accumulates run outcomes
for AI-driven parameter suggestions and walltime estimation.

At the start of any run-related conversation:

1. Read `run_journal.jsonl` if it exists
2. Summarize relevant past runs for the current discussion
3. When proposing parameters, cite relevant history:
   "Based on run 20260217 (8x8 D=8 PBC), SR lr=0.1 converged in 22
   iters with 1 spike. Suggesting lr=0.05."

After fetching results (`fetch` command auto-appends to journal):

4. Review the new journal entry
5. Analyze convergence quality (see checklist below)
6. Suggest improvements if warranted

## Post-Fetch Analysis Checklist

After fetching results, update `run_history.json` first:
1. Set job `status` to COMPLETED or FAILED
2. Fill in `ended`, `last_energy`, `last_energy_per_site`
3. If failed, write `failure_reason` (diagnosis, not just "killed")
4. Update `current_status` and `best_energy_per_site`

Then evaluate:

- **Energy trend**: converging / plateaued / diverging?
- **Spike count**: 0 ideal, >2 suggests lr or CGDiagShift issue
- **Accept rate**: 0.2-0.3 is healthy for spin-1/2
- **Walltime utilization**: <10% means walltime too generous, suggest 5x median
- **CG iterations**: hitting CGMaxIter often? suggest increasing it
- **Gradient norm**: should decrease; flat means optimization stuck
- **EvalT / UpdateT**: compare to similar past runs, flag regressions

## Improvement Suggestion Categories

When patterns emerge across multiple runs, suggest:

**Parameter tuning:**
- "D=8 runs consistently need lr<=0.05 to avoid spikes"
- "WarmUp=128 insufficient for 12x12; energy unstable in first iterations"

**HPC resource optimization:**
- "256G56c uses only 30GB for D=5; smaller partition available?"
- "Last 3 runs used <7% walltime; suggest 5x median = 2d12h"
- "Consider --ntasks 28 with OMP_NUM_THREADS=2 for memory-bound D=10"

**Code improvements:**
- "SpikeFactorErr=100 triggers in 4/5 D>=8 runs; consider tightening default
  in enhanced_params_parser.h"
- "CGDiagShift=1e-4 causes CG stalling at D=10; scale with D?"

**Performance:**
- "EvalT increased 3x from D=5 to D=8 but UpdateT stayed flat;
  contraction is the bottleneck, not the CG solver"

## Run History Protocol

Every cluster run date-subdirectory must have a `run_history.json`. This is
the single source of truth for what happened, why, and the current state.

### Schema

```json
{
  "purpose": "One-line goal of this run",
  "context": "Why this approach was chosen, what alternatives were tried",
  "pipeline": ["su", "tile", "vmc", "measure"],
  "physics": {
    "model": "SquareHeisenberg",
    "Lx": 8, "Ly": 8, "J2": 0.0, "D": 5, "BC": "Periodic"
  },
  "source_data": "relative path to initial tpsfinal, or null if from SU",
  "jobs": [
    {
      "id": 522405,
      "stage": "vmc",
      "status": "COMPLETED|FAILED|RUNNING|PENDING|CANCELLED",
      "submitted": "ISO 8601",
      "ended": "ISO 8601 or null",
      "failure_reason": "human-readable diagnosis or null",
      "iterations_completed": 30,
      "last_energy": -42.999,
      "last_energy_per_site": -0.672,
      "params_changed": {"InitialStepSelectorEnabled": false},
      "notes": "free-text context"
    }
  ],
  "current_status": "completed|running|failed|staged",
  "best_energy_per_site": -0.672,
  "best_energy_error": 0.00003
}
```

### Update rules

- **On submission**: Add job entry with `status: "RUNNING"`.
- **On completion**: Update `status`, `ended`, energy results.
- **On failure**: Set `status: "FAILED"`, write `failure_reason` with
  root cause diagnosis (not just "killed" or "crashed").
- **On resubmission**: Record `params_changed` in the new job entry.
- **Prune transient retries**: If a job was just a stepping stone to
  finding the real fix, delete its log files and don't record it.
  Only keep meaningful attempts that contributed to diagnosis or results.

### Discovery

To scan all runs without reading logs:

```bash
# Find all documented runs
glob("run/**/run_history.json")

# Quick status overview
jq '.current_status' run/**/run_history.json

# Find all completed runs with energy
jq 'select(.current_status=="completed") | .best_energy_per_site' run/**/run_history.json
```

### Failure signatures to document

When writing `failure_reason`, include actionable diagnosis:
- "InitialStepSelector bumped LR from 0.1 to 0.3, causing CG divergence"
- "Killed by cluster disk failure and reboot on 2026-02-21"
- "Loop update unstable: E0 drifted upward, advanced stop never triggered"

Not: "job failed", "crashed", "exit code 127"

## Benchmark Energies

Reference energies for validating results are stored in
`data/benchmarks/qmc_square_heisenberg.json`. When plotting energy
trajectories, read this file and overlay the appropriate QMC benchmark
line for the matching model, J2, and lattice size.

## Walltime Strategy

- **No history**: Default 60 days (safe ceiling for unknown jobs)
- **With history**: AI reads journal, finds similar runs (same lattice/D/BC),
  suggests 5x median walltime. Displayed in summary.
- **User override**: `--walltime` always wins
