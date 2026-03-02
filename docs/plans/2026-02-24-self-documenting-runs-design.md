# Design: Self-Documenting Run Directories

**Date**: 2026-02-24
**Status**: Approved

## Problem

Cluster run directories accumulate history (failed jobs, parameter changes,
restarts) that only the person who managed them understands. Future AI agents
or collaborators looking at the directories must grep through multiple slurm
logs to reconstruct what happened, why, and what the current state is.

## Solution

Add a machine-readable `run_history.json` to each date-subdirectory under
`run/`. This file is the single source of truth for that run's purpose,
physics, job history, and current status.

## Schema

```json
{
  "purpose": "Brief description of what this run aims to accomplish",
  "context": "Strategic context: why this approach, what alternatives were tried",
  "pipeline": ["vmc", "measure"],
  "physics": {
    "model": "SquareHeisenberg",
    "Lx": 8, "Ly": 8, "J2": 0.0, "D": 5, "BC": "Periodic"
  },
  "source_data": "relative path to initial tpsfinal, or null if from SU",
  "jobs": [
    {
      "id": 522224,
      "stage": "vmc",
      "status": "COMPLETED|FAILED|RUNNING|PENDING|CANCELLED",
      "submitted": "2026-02-23T11:39:00",
      "ended": "2026-02-23T14:03:00",
      "failure_reason": "InitialStepSelector bumped LR from 0.1 to 0.3, causing CG divergence",
      "iterations_completed": 1,
      "last_energy": -42.999,
      "last_energy_per_site": -0.672,
      "params_changed": {"InitialStepSelectorEnabled": false},
      "notes": "Restored clean Yubin tpsfinal after this failure"
    }
  ],
  "current_status": "completed|running|failed|staged",
  "best_energy_per_site": -0.672,
  "best_energy_error": 0.00003
}
```

### Field descriptions

- **purpose**: One-line description of the run's goal.
- **context**: Why this approach was chosen. E.g., "Yubin tensor gives better
  E/site than any SU/loop-update initial state."
- **pipeline**: Ordered list of stages this run executes.
- **physics**: Key physical parameters extracted from physics JSON.
- **source_data**: Where the initial TPS came from (relative to `run/`).
- **jobs**: Array of meaningful job attempts. Only record attempts that
  contributed to diagnosis or produced results. Delete transient retries
  and their logs.
- **current_status**: High-level summary for quick scanning.
- **best_energy_per_site / best_energy_error**: Final or best result.

### Jobs array conventions

- Append-only for meaningful attempts (not every transient retry).
- **params_changed**: Documents what was tuned vs. the previous attempt.
- **failure_reason**: Human-readable diagnosis so future agents don't
  need to grep logs.
- Transient retries (e.g., resubmission before finding root cause) should
  be deleted along with their log files.

## Update discipline

- **On submission**: Add job entry with `status: "RUNNING"`.
- **On completion**: Update status, energy results.
- **On failure + diagnosis**: Update `failure_reason`, then record the fix
  in the next job entry's `params_changed`.
- **On pruning**: Remove transient job entries and their cluster log files.

## Discovery

A new agent arriving at the repo:
1. `glob("run/**/run_history.json")` — finds all documented runs.
2. Read any file for full context without touching logs.
3. `jq '.current_status' run/**/run_history.json` — quick status scan.

## Scope

Applies to all cluster run date-subdirectories. Local test runs (the
`20260214_*` dirs) already have `status.txt` and `run_info.txt` and are
out of scope.

## Future extensions

- Top-level `run/index.json` for fast overview (generate from individual files).
- README auto-generation from JSON.
- Integration with `cluster_run.py` to auto-create/update `run_history.json`.
