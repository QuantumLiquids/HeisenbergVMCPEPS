# Changelog: Upstream PEPS API Changes

## PEPS main 68452fe - Row/column axis update - Applied 2026-10-02

Upstream PEPS `main` adds `qlpeps::MCUpdateSquareAxisAutoregressiveOBC`, a rejection-free
row/column ("axis") autoregressive sweep updater for finite OBC states, the composite
`qlpeps::MCUpdateSequence`, and `qlpeps::LabelCountTable`, which keeps a count that dense
tensors do not carry (here N_up) fixed during sampling. `vmc_optimize` and `mc_measure` can now
run it before the local `MCUpdateSquareTNN3SiteExchangeOBC` sweep.

### Why it matters

A row or column move redraws a whole slice from its boundary-MPS window, which can shorten
autocorrelation times where the local 3-site exchange mixes slowly. It is not cheap. Besides
sampling every slice, the axis pass performs its own BMPS compressions, and its row pass
consumes the lower boundary MPS as it moves down, which the 3-site sweep then regrows. On the 12x12
D=8 Heisenberg fixtures with `Dbmps_max = MCAxisUpdateDmax = 32`, a sweep with the axis update
took 2.9 times (dense with the count table) and 3.8 times (`-DU1SYM`) as long as a 3-site sweep
alone, with 77 instead of 33 BMPS compressions (PEPS `profiler/README.md`, section "Measured",
workload `heisenberg-sector`; per-sample figures in `tutorials/04-parameter-reference.md`,
section 4.8). Compare autocorrelation per wall-clock time.

### New JSON keys (VMC and measure algorithm files)

- `MCAxisUpdate` — optional JSON bool, default `false`. When `true`, every OBC sweep runs
  `MCUpdateSequence<MCUpdateSquareAxisAutoregressiveOBC, MCUpdateSquareTNN3SiteExchangeOBC>`.
- `MCAxisUpdateDmin` (default `1`), `MCAxisUpdateDmax` (default `Dbmps_max`),
  `MCAxisUpdateTruncErr` (default `0.0`) — the slice-MPS truncation
  (`qlpeps::BMPSSliceSamplerParams`). No Metropolis-Hastings key: the update is rejection-free,
  so scan `MCAxisUpdateDmax` upward.

With `MCAxisUpdate=true` the parsers throw `std::invalid_argument` for PBC, nonzero
`SpinInversionParity`, `MCAxisUpdateDmax == 0`, `MCAxisUpdateDmin` outside
`[1, MCAxisUpdateDmax]`, `MCAxisUpdateTruncErr` outside `[0, 1)`, and, when S_z is conserved
(`-DU1SYM` or `MCRestrictU1=true`), `MCAxisUpdateDmax < floor(max(Lx, Ly) / 2) + 1` (PEPS
design section 5.4). On dense tensors `MCRestrictU1=true` adds the N_up table `{{1},{0}}`;
the axis update is the only sampler that honours `MCRestrictU1`. Logs, the dense versus
`-DU1SYM` table and the hidden-S_z caveat: `tutorials/04-parameter-reference.md`, section 4.8.

### Build

Requires PEPS `main` at or after `68452fe` (the row/column axis update series). PEPS still
reports version 0.2.2 there, as before these headers existed, so
`find_package(PEPS 0.2.2 CONFIG REQUIRED)` also accepts an older 0.2.2 install;
`src/model_updater_factory.h` then stops the build with an `#error` that names the
requirement. The version requirement moves up with the next PEPS release.

### Backward compatibility

Off by default. With `MCAxisUpdate` absent or `false`, OBC runs use the same
default-constructed `MCUpdateSquareTNN3SiteExchangeOBC` and print the same logs as before;
PBC and spin-inversion runs are unchanged. Existing parameter files need no change.
`RunMeasureByModel` gained `(axis_update, mc_params)` arguments. New fast tests:
`test_axis_update_params` (keys, validation, mapping, logs) and `test_axis_update_sweep`
(sweeps of the enabled updater on a small simple-update state of this build's types).

---

## PEPS v0.2.1 - HOTRG PBC contractor - Applied 2026-09-21

Upstream PEPS 0.2.1 adds `qlpeps::HOTRGContractor`, a second periodic-boundary
contraction backend alongside `qlpeps::TRGContractor`. It plugs into
`VMCPEPSOptimizer` and `MCPEPSMeasurer` at the same template slot, and the
existing PBC updater and J1-J2 XXZ PBC solver accept either one.

### Why it matters

TRG only accepts a square lattice with linear size `2^k` or `3*2^k`. HOTRG
accepts any `Ly x Lx` with both dimensions at least 2, so rectangles and odd
sizes such as 5x5 are now reachable under PBC. Accuracy guidance is unchanged:
the environment bond dimension should be at least `D^2`.

### New JSON keys (VMC and measure algorithm files)

- `PBCContractor` — optional, `TRG` (default) or `HOTRG`, case-insensitive.
- `HOTRGDmin`, `HOTRGDmax`, `HOTRGTruncErr` — required when `PBCContractor` is
  `HOTRG`. HOTRG inverts nothing, so there is no counterpart to
  `TRGInvRelativeEps`.

Selecting `HOTRG` without those three keys throws
`PBC with PBCContractor=HOTRG requested but HOTRG params are missing in
algorithm JSON. Require: HOTRGDmin, HOTRGDmax, HOTRGTruncErr.`
Any other `PBCContractor` value throws `PBCContractor must be TRG or HOTRG.`

### Build

`CMakeLists.txt` now requires `find_package(PEPS 0.2.1 CONFIG REQUIRED)`.
Point `-DPEPS_DIR` at a PEPS 0.2.1 install.

### Backward compatibility

`PBCContractor` defaults to `TRG`, so existing PBC parameter files keep
selecting TRG and produce the same runs as before. No OBC path changed.

### Example parameter files

- `params/quickstart/vmc_local_2x2_pbc_hotrg_n1.json`
- `params/quickstart/measure_local_2x2_pbc_hotrg_n1.json`

---

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
