# Changelog: Upstream PEPS API Changes

## PEPS commit e572888 - Applied 2026-03-03

Cluster binary rebuilt 2026-03-02 20:58. **Jobs submitted before this date
do NOT have these features.**

### JSONL structured logging (auto-enabled)

Per-iteration machine-readable log auto-generated at
`vmc/energy/optimization_log.jsonl`. No configuration needed — upstream
auto-configures when `tps_dump_base_name` is set (always true in our runs).

### PeriodicStepSelector fixes

- Selector **no longer triggers at iter 0** (previously caused 2x EvalT
  overhead on the very first iteration).
- New `SelectorT` field in log output separates selector time from UpdateT.
- Log line now shows `CG resid` field.

### Affected jobs

Jobs using the old binary (before 2026-03-02 rebuild):
- 523623, 523634, 524096, 524117, 524119, 524679 — no JSONL, no SelectorT,
  iter 0 has 2x EvalT spike in UpdateT.

Jobs using the new binary (after rebuild): will have JSONL and fixes.

---

## PEPS post-v0.1.0 (up to f13022c) - Applied 2026-02-28

### ConjugateGradientParams (aggregate, no constructors)

Old: `ConjugateGradientParams(max_iter, tol, restart, diag_shift)`
New: designated-init aggregate with fields:
- `.max_iter` (size_t)
- `.relative_tolerance` (double) — was `.tolerance`, now in norm-space (old was squared-residual)
- `.absolute_tolerance` (double, default 0.0)
- `.residual_recompute_interval` (int) — was `.residue_restart_step`
- `.orthogonality_threshold` (double, default 0.5)

`diag_shift` removed from CG; moved to `StochasticReconfigurationParams`.

### StochasticReconfigurationParams (aggregate, no constructors)

Old: `StochasticReconfigurationParams(cg_params, normalize_update)`
New: designated-init aggregate with fields:
- `.cg_params` (ConjugateGradientParams)
- `.diag_shift` (double, default 0.001) — moved from CG params
- `.normalize_update` (bool)
- `.adaptive_diagonal_shift` (double, default 0.0)

### Step selector rename

- `AutoStepSelectorParams` -> `PeriodicStepSelectorParams`
- `BaseParams::auto_step_selector` -> `BaseParams::periodic_step_selector`
- Builder: `.SetAutoStepSelector(...)` -> `.SetPeriodicStepSelector(...)`

### JSON parameter key aliases

Old JSON keys are still accepted with deprecation warnings. New keys align with upstream names.

| Old key (deprecated) | New key | Conversion |
|---|---|---|
| `CGTol` | `CGRelativeTolerance` | auto `sqrt()` (old was squared-residual, new is norm-space) |
| `CGResidueRestart` | `CGResidualRecomputeInterval` | name only |
| `CGDiagShift` | `SRDiagShift` | moved from CG to SR scope |

Both old and new keys are present → error (ambiguous).

### MinSR optimizer (commit 6546955) - Wired 2026-03-01

New optimizer variant: MinSR (Minimum-step Stochastic Reconfiguration, Chen & Heyl 2024).
Solves Ns x Ns system instead of Np x Np when Np >> Ns.

`MinSRParams` struct (has constructors):
- `.r_pinv` (double, default 1e-12) — relative pseudo-inverse cutoff
- `.a_pinv` (double, default 0.0) — absolute pseudo-inverse cutoff
- `.soft_cutoff` (bool, default true) — soft cutoff formula (Eq. 23)
- `.solver_mode` (MinSRSolverMode, default kAuto) — Auto/Replicated/Distributed

`MinSRSolverMode` enum: `kAuto`, `kReplicated` (LAPACK), `kDistributed` (ScaLAPACK)

JSON keys in HeisenbergVMCPEPS:
- `"OptimizerType": "MinSR"`
- `MinSRRPinv` (default 1e-12)
- `MinSRAPinv` (default 0.0)
- `MinSRSoftCutoff` (default true)
- `MinSRSolverMode` (default "Auto"; accepts Auto/Replicated/Distributed)

### Other changes (not affecting this codebase directly)

- `step_length_trajectory` -> `learning_rate_trajectory`
- `BoundedGradientUpdate` removed
- `CGResult.converged` field -> `CGTerminationReason reason` + `converged()` method
- `SRSMatrix` ctor: 3 args -> 4 args (added `MPI_Comm`)
- `InplaceMultiplyMPO()` -> `MultiplyMPOInplace()`
