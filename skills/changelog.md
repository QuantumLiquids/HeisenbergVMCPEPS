# Changelog: Upstream PEPS API Changes

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

### Other changes (not affecting this codebase directly)

- `step_length_trajectory` -> `learning_rate_trajectory`
- `BoundedGradientUpdate` removed
- `CGResult.converged` field -> `CGTerminationReason reason` + `converged()` method
- `SRSMatrix` ctor: 3 args -> 4 args (added `MPI_Comm`)
- `InplaceMultiplyMPO()` -> `MultiplyMPOInplace()`
