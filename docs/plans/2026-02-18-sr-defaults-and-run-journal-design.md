# SR Default Profiles & Run Journal Learning System

Date: 2026-02-18
Status: Approved
Builds on: `docs/plans/2026-02-17-cluster-workflow-skill-cli-design.md`

## 1. Overview

Two improvements to the cluster workflow CLI and skill:

1. **SR parameter defaults**: Battle-tested universal CG constants with auto-scaled
   WarmUp. No named profiles — just smart defaults displayed to the user.
2. **Run journal learning system**: A persistent `run_journal.jsonl` that accumulates
   run outcomes, enabling the AI to propose better parameters, walltime estimates,
   and code improvements based on history.

## 2. SR Parameter Defaults

### 2.1 Universal Constants

These SR/CG parameters are constant across all system sizes:

| Parameter              | Value | Rationale                                    |
|------------------------|-------|----------------------------------------------|
| CGMaxIter              | 200   | Upper bound; CG usually converges in 10-50   |
| CGTol                  | 1e-6  | Standard convergence tolerance                |
| CGResidueRestart       | 10    | Restart CG every 10 iters (20 for large)     |
| CGDiagShift            | 1e-4  | Fisher matrix regularization                 |
| NormalizeUpdate        | false | Standard SR practice                         |
| MCLocalUpdateSweeps    | 1     | Single sweep between samples                 |
| LearningRate           | 0.1   | Default; overridable via `--lr`              |

CGResidueRestart increases to 20 when `max(Lx, Ly) >= 16`.

### 2.2 Auto-Scaled WarmUp

```
WarmUp = max(100, min(500, Lx * Ly * 2))
```

Examples:
- 4x4 (16 sites): WarmUp = 100 (clamped to min)
- 8x8 (64 sites): WarmUp = 128
- 12x12 (144 sites): WarmUp = 288
- 16x16 (256 sites): WarmUp = 500 (clamped to max)

Displayed in `new-case` summary as:
```
WarmUp:        128 (auto: 2 x Lx x Ly, range [100, 500])
```

Overridable via `--override WarmUp=300`.

### 2.3 QoS Support

Add `#SBATCH --qos=wanghxqos` to all slurm templates. CLI flag `--qos` with
default `wanghxqos`. MaxJobs=30 allows 10 simultaneous 3-stage pipelines.

### 2.4 Walltime Strategy

Three tiers:

1. **No history (default)**: `60-00:00:00` (2 months). Safe ceiling for first runs
   of unknown duration.
2. **With history**: AI reads journal, finds similar past runs (same lattice/D/BC),
   computes `5 x median_walltime_used`, and suggests that as the limit. Displayed:
   `Walltime: 2-12:00:00 (auto: 5x median of 3 similar runs)`.
3. **User override**: `--walltime` always wins.

## 3. Run Journal

### 3.1 File Location

`run_journal.jsonl` in repo root. Gitignored — this is local operational state,
not committed science data.

### 3.2 Entry Schema

Each line is a JSON object appended after `fetch`:

```json
{
  "case": "8x8J2=0D8",
  "run_id": "20260217_211901",
  "timestamp": "2026-02-17 21:19:01",
  "params": {
    "lx": 8, "ly": 8, "D": 8, "J2": 0.0, "bc": "pbc",
    "optimizer": "SR", "lr": 0.1, "iterations": 30,
    "vmc_samples": 11200, "warmup": 200,
    "cg_max_iter": 200, "cg_tol": 1e-6, "cg_diag_shift": 1e-4
  },
  "outcome": {
    "converged": true,
    "final_energy": -25.134,
    "final_energy_err": 0.012,
    "best_energy": -25.201,
    "best_iter": 22,
    "gradient_norm_final": 0.08,
    "spike_count": 1,
    "accept_rate": 0.24,
    "total_iterations": 30
  },
  "hpc": {
    "partition": "256G56c",
    "qos": "wanghxqos",
    "ntasks": 56,
    "walltime_requested": "60-00:00:00",
    "walltime_used_seconds": 18432,
    "utilization": 0.004
  },
  "performance": {
    "eval_time_mean": 15.6,
    "update_time_mean": 1.9,
    "total_time_per_iter": 17.5
  },
  "insights": []
}
```

### 3.3 Lifecycle

```
run completes
    |
    v
fetch command --> parse energy_trajectory.csv + slurm logs
    |
    v
append structured entry to run_journal.jsonl
    |
    v
AI analyzes outcome --> adds insights to entry
    |
    v
new run requested --> AI reads journal
    |                  consults similar past runs
    |                  proposes informed parameters
    v
AI notices cross-run patterns --> suggests improvements
    (params / HPC / C++ code / performance)
```

### 3.4 Journal Parsing (in fetch command)

The `fetch` command gains a post-processing step:

1. Read `energy/energy_trajectory.csv` columns: iteration, energy, energy_err,
   gradient_norm
2. Compute: final_energy, best_energy, best_iter, gradient_norm_final,
   converged (energy_err < 0.1 * |energy| in last 5 iters)
3. Count spike recoveries from `energy/spike_trigger.csv` if present
4. Parse slurm log for: walltime used, accept rate
5. Parse slurm log for: EvalT, UpdateT per iteration (performance metrics)
6. Build the journal entry and append to `run_journal.jsonl`

### 3.5 What the AI Does NOT Do Automatically

- Does not autonomously monitor jobs or fetch results
- Does not modify parameters without user approval
- Does not commit journal entries (the file is gitignored)
- Does not push code changes — only suggests them in conversation

## 4. Skill Enhancements

New sections added to `skills/vmcpeps-cluster-workflow.md`:

### 4.1 Run Journal Protocol

```
At the start of any run-related conversation:
1. Read run_journal.jsonl if it exists
2. Summarize relevant past runs for the current discussion
3. When proposing parameters, cite relevant history:
   "Based on run 20260217 (8x8 D=8 PBC), SR lr=0.1
    converged in 22 iters with 1 spike. Suggesting lr=0.05."
```

### 4.2 Post-Fetch Analysis Checklist

```
After fetching results, evaluate:
- Energy trend: converging / plateaued / diverging?
- Spike count: 0 ideal, >2 suggests lr or CGDiagShift issue
- Accept rate: 0.2-0.3 is healthy for spin-1/2
- Walltime utilization: <10% means walltime too generous
- CG iterations per step: hitting CGMaxIter often? increase it
- Gradient norm trend: should decrease; flat = optimization stuck
- EvalT/UpdateT: compare to similar past runs, flag regressions
```

### 4.3 Improvement Suggestion Categories

```
When patterns emerge across multiple runs, suggest:

Parameter tuning:
  "D=8 runs consistently need lr<=0.05 to avoid spikes"

HPC resource optimization:
  "256G56c uses only 30GB for D=5; smaller partition available?"
  "Last 3 runs used <7% walltime; suggest 5x median = 2d12h"

Code improvements:
  "SpikeFactorErr=100 triggers in 4/5 D>=8 runs — consider PR
   to change default in enhanced_params_parser.h"
  "CGDiagShift=1e-4 causes CG stalling at D=10; scale with D?"

Performance:
  "EvalT increased 3x from D=5 to D=8 but UpdateT stayed flat —
   contraction is the bottleneck, not the CG solver"
  "OMP_NUM_THREADS=1 with 56 ranks; could try 2 threads x 28 ranks"
```

## 5. CLI Changes Summary

### 5.1 cluster_run.py Modifications

- `build_params()`: Replace hardcoded WarmUp=100 with auto-scaling formula.
  Add CGResidueRestart scaling for large systems.
- `generate_slurm_*()`: Add `#SBATCH --qos={qos}` line. Change default
  walltime from `72:00:00` to `60-00:00:00`.
- `print_summary()`: Show WarmUp with "(auto: ...)" hint. Show walltime
  with journal-based suggestion if available.
- `cmd_fetch()`: After downloading results, call `append_journal_entry()`
  to parse and append to `run_journal.jsonl`.
- New function: `append_journal_entry(case_dir, meta)` — parses CSVs and
  logs, builds journal entry, appends to jsonl.
- New function: `read_journal()` — loads and returns all journal entries.
- New function: `suggest_walltime(journal, lx, ly, d, bc)` — finds similar
  runs and returns 5x median walltime.
- New CLI flag: `--qos` (default: `wanghxqos`)
- Changed default: `--walltime` from `72:00:00` to `60-00:00:00`

### 5.2 .gitignore Addition

```
run_journal.jsonl
```

### 5.3 Test Additions

- `test_compute_warmup()`: verify scaling formula and clamping
- `test_cg_residue_restart_scaling()`: 10 for small, 20 for large
- `test_append_journal_entry()`: mock CSV/log parsing, verify schema
- `test_suggest_walltime()`: with/without history, verify 5x median
- `test_qos_in_slurm_output()`: verify `#SBATCH --qos` appears
