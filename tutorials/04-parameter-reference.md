## Parameter Reference (by Command)

Use this page as the key-by-key contract.

For runnable minimal examples, see `tutorials/03-recipes.md`.

### 1) Shared Physics File: `physics_params.json`

Required keys:

- `Lx` (int)
- `Ly` (int)
- `J2` (double)
- `ModelType` (string): `SquareHeisenberg`, `SquareXY`, `TriangleHeisenberg`

Optional keys:

- `BoundaryCondition` (string): `Open`/`OBC`/`Periodic`/`PBC` (case-insensitive)
  - default: `Open`
- `RemoveCorner` (bool, legacy compatibility key)
  - ignored by unified square/triangle drivers

Runtime effects:

- `ModelType` selects solver dispatch path.
- `BoundaryCondition` selects contraction backend: OBC -> BMPS, PBC -> TRG.

### 2) `simple_update_algorithm_params.json`

Required keys:

- `Tau` (double)
- `Step` (int)
- `Dmin` (int)
- `Dmax` (int)
- `TruncErr` (double)
- `ThreadNum` (int)

Optional advanced-stop keys:

- `AdvancedStopEnabled` (bool, default `false`)
- `AdvancedStopEnergyAbsTol` (double, default `1e-8`)
- `AdvancedStopEnergyRelTol` (double, default `1e-10`)
- `AdvancedStopLambdaRelTol` (double, default `1e-6`)
- `AdvancedStopPatience` (int, default `3`)
- `AdvancedStopMinSteps` (int, default `10`)

Optional tau-schedule keys:

- `TauScheduleEnabled` (bool, default `false`)
- `TauScheduleTaus` (string, required when `TauScheduleEnabled=true`)
  - comma-separated positive doubles (example: `"0.5,0.2,0.1,0.05,0.02"`)
- `TauScheduleStepCaps` (string, optional)
  - comma-separated positive integers
  - if omitted, every stage uses global `Step`
- `TauScheduleRequireConverged` (bool, default `true`)
  - when true, a stage is considered failed if `GetLastRunSummary().converged=false`
  - requires advanced-stop to be active (`AdvancedStopEnabled=true` or any `AdvancedStop*` tuning keys present)
- `TauScheduleDumpEachStage` (bool, default `false`)
- `TauScheduleDumpDir` (string, default `"tau_schedule"`)
- `TauScheduleAbortOnStageFailure` (bool, default `true`)

Advanced-stop activation:

- Enabled when `AdvancedStopEnabled=true`.
- Also auto-enabled when any advanced-stop tuning key is present.
- If `AdvancedStopEnabled=false` is set explicitly, advanced-stop is disabled even if tuning keys are present.

Advanced-stop convergence rule:

- Gate condition is `energy AND lambda`, then apply `patience` and `min_steps`.
- Energy criterion:
  - `|ΔE| <= max(energy_abs_tol, energy_rel_tol * max(1, |E_prev|, |E_curr|))`
- Lambda criterion:
  - compute per-bond relative L2 drift on lambda diagonals and take the global maximum;
  - require `max_lambda_drift <= lambda_rel_tol`.
- If bond dimensions change between two sweeps, lambda drift is skipped and the convergence streak resets.
- Stop when gate passes for `AdvancedStopPatience` consecutive sweeps and executed sweeps >= `AdvancedStopMinSteps`.

Runtime effects:

- `Tau` and `Step` control imaginary-time evolution length/granularity.
- `Dmin`/`Dmax` set SU bond-dimension range.
- Advanced-stop (when active) can terminate before `Step`; otherwise run always executes to `Step`.
- Output always targets `tpsfinal/` (SITPS) and `peps/`.
- When `TauScheduleEnabled=false` (default), driver runs one stage using global `Tau` + `Step` (legacy behavior).
- When `TauScheduleEnabled=true`, stage `tau`/`step_cap` override global `Tau`/`Step` at runtime.
- Tau stages run in listed order on the same evolving PEPS.
- Driver writes machine-readable schedule summaries to:
  - `<TauScheduleDumpDir>/schedule_summary.json`
  - `<TauScheduleDumpDir>/schedule_summary.csv`
- If `TauScheduleDumpEachStage=true`, driver also dumps stage snapshots under:
  - `<TauScheduleDumpDir>/stage_XX_tau_<value>/tpsfinal`
  - `<TauScheduleDumpDir>/stage_XX_tau_<value>/peps`

### 3) `loop_update_algorithm_params.json`

Required keys:

- `Tau` (double)
- `Step` (int)
- `Dmin` (int)
- `Dmax` (int)
- `TruncErr` (double)
- `ThreadNum` (int)

Optional advanced-stop keys (same semantics as simple update):

- `AdvancedStopEnabled` (bool, default `false`)
- `AdvancedStopEnergyAbsTol` (double, default `1e-8`)
- `AdvancedStopEnergyRelTol` (double, default `1e-10`)
- `AdvancedStopLambdaRelTol` (double, default `1e-6`)
- `AdvancedStopPatience` (int, default `3`)
- `AdvancedStopMinSteps` (int, default `10`)

Optional loop-truncation keys:

- `ArnoldiTol` (double, default `1e-10`)
- `ArnoldiMaxIter` (int, default `200`)
- `LoopInvTol` (double, default `1e-8`)
- `FETTolerance` (double, default `1e-12`)
- `FETMaxIter` (int, default `30`)
- `CGMaxIter` (int, default `100`)
- `CGTol` (double, default `1e-10`)
- `CGResidueRestart` (int, default `20`)
- `CGDiagShift` (double, default `0.0`)

Current scope constraints (hard-fail if violated):

- `ModelType` must be `SquareHeisenberg`
- `J2` must be exactly zero within numerical tolerance

Runtime effects:

- Driver runs loop-update sweeps with `qlpeps::LoopUpdateParams` built from this JSON.
- Advanced-stop (when active) can terminate before `Step`; otherwise run executes to `Step`.
- Output always targets `tpsfinal/` (SITPS) and `peps/`.
- Driver logs include an optional advanced-stop summary:
  - `converged`, `stop_reason`, `executed_steps/Step`.

### 4) `vmc_algorithm_params.json`

#### 4.1 Required baseline keys

- `MC_total_samples` (int)
- `WarmUp` (int)
- `MCLocalUpdateSweepsBetweenSample` (int)
- `ThreadNum` (int, optional in parser, default `1`, but should be set explicitly)
- `OptimizerType` (optional in parser, default `StochasticReconfiguration`)
- `MaxIterations` (optional in parser, default `10`)
- `LearningRate` (optional in parser, default `0.01`)

Important SR note:

- If `OptimizerType` is SR / `StochasticReconfiguration`, these become required:
  - `CGMaxIter`
  - `CGTol`
  - `CGResidueRestart`
- `CGDiagShift` default is `0.0`.
- `NormalizeUpdate` default is `false`.

#### 4.2 Backend keys (selected by boundary condition)

OBC/BMPS:

- required: `Dbmps_max`
- optional:
  - `Dbmps_min` (default `Dbmps_max`)
  - `BMPSTruncErr` (default `0.0`)
  - `MPSCompressScheme` (default `SVD`)
  - `BMPSConvergenceTol` (used by variational schemes)
  - `BMPSIterMax` (used by variational schemes)

PBC (backend selected by `PBCContractor`):

- `PBCContractor` (optional, default `TRG`; accepts `TRG` or `HOTRG`,
  case-insensitive)

PBC/TRG (`PBCContractor` absent or `TRG`):

- required: `TRGDmin`, `TRGDmax`, `TRGTruncErr`
- optional: `TRGInvRelativeEps` (default `1e-12`)
- geometry: square lattice with linear size `2^k` or `3*2^k`

PBC/HOTRG (`"PBCContractor": "HOTRG"`):

- required: `HOTRGDmin`, `HOTRGDmax`, `HOTRGTruncErr`
- no inversion parameter: HOTRG inverts nothing, so `TRGInvRelativeEps` has no
  HOTRG counterpart
- geometry: any `Ly x Lx` with both dimensions at least 2

#### 4.3 IO keys

- `WavefunctionBase` (string, default `"tps"`)
  - load path is `WavefunctionBase + "final"` -> usually `tpsfinal/`
- `ConfigurationLoadDir` (string, default `WavefunctionBase + "final"`)
- `ConfigurationDumpDir` (string, default `WavefunctionBase + "final"`)

Runtime effects:

- `configuration{rank}` is loaded from `ConfigurationLoadDir` when available.
- final configuration is dumped to `ConfigurationDumpDir`.

#### 4.4 Optimizer-specific keys

SGD:

- `Momentum` (default `0.0`)
- `Nesterov` (default `false`)
- `WeightDecay` (default `0.0`)

Adam:

- `Beta1` (default `0.9`)
- `Beta2` (default `0.999`)
- `Epsilon` (default `1e-8`)
- `WeightDecay` (default `0.0`)

AdaGrad:

- `Epsilon` (default `1e-8`)
- `InitialAccumulator` (default `0.0`)

LBFGS:

- `LBFGSHistorySize` (default `10`)
- `LBFGSToleranceGrad` (default `1e-5`)
- `LBFGSToleranceChange` (default `1e-9`)
- `LBFGSMaxEval` (default `20`)
- `LBFGSStepMode` (default `Fixed`; accepts `Fixed`, `StrongWolfe`, `kFixed`, `kStrongWolfe`)
- `LBFGSWolfeC1` (default `1e-4`)
- `LBFGSWolfeC2` (default `0.9`)
- `LBFGSMinStep` (default `1e-8`)
- `LBFGSMaxStep` (default `1.0`)
- `LBFGSMinCurvature` (default `1e-12`)
- `LBFGSUseDamping` (default `true`)
- `LBFGSMaxDirectionNorm` (default `1e3`)
- `LBFGSAllowFallbackToFixedStep` (default `false`)
- `LBFGSFallbackFixedStepScale` (default `0.2`)

#### 4.5 Step selectors (SGD/SR only)

- `InitialStepSelectorEnabled` (default `false`)
- `InitialStepSelectorMaxLineSearchSteps` (default `3`)
- `InitialStepSelectorEnableInDeterministic` (default `false`)
- `AutoStepSelectorEnabled` (default `false`)
- `AutoStepSelectorEveryNSteps` (default `10`)
- `AutoStepSelectorPhaseSwitchRatio` (default `0.3`)
- `AutoStepSelectorEnableInDeterministic` (default `false`)

Constraints:

- Selectors only valid for `SGD` or SR.
- Selectors cannot be combined with `LRScheduler`.
- `LearningRate` must be positive when any selector is enabled.
- `InitialStepSelectorMaxLineSearchSteps > 0` when initial selector enabled.
- `AutoStepSelectorEveryNSteps > 0` and `AutoStepSelectorPhaseSwitchRatio in [0,1]` when auto selector enabled.

#### 4.6 Misc optional knobs

- `LRScheduler`: `ExponentialDecay`, `CosineAnnealing`, `Plateau`
- `ClipNorm`, `ClipValue`
- spike recovery keys (`Spike*`)
- checkpoint keys (`CheckpointEveryNSteps`, `CheckpointBasePath`)
- `MCRestrictU1` (default `true`). The local updaters always conserve S_z, so the key
  matters only with the opt-in axis update (section 4.8), where it selects the S_z count
  table for dense tensors.
- `InitialConfigStrategy` (default `Random`; accepts `Random`, `Neel`, `ThreeSublatticePolarizedSeed`)

#### 4.7 Spin-inversion projection (opt-in)

`SpinInversionParity` is an integer in the VMC algorithm file: `0` (default)
keeps the plain PEPS path; `+1` or `-1` optimizes the coherent amplitude
`psi(x) + parity * psi(Fx)`, where `F` swaps up/down at every site.
This is spin inversion, not a singlet projection. Parameters are shared between
both branches; sampling, energies, gradients, SR and MinSR use their coherent sum.

The initial implementation requires `SquareHeisenberg` or `SquareXY`, OBC,
`J2=0`, no removed corners, `MCRestrictU1=true`, `Lx,Ly >= 2`, and even `Lx*Ly`.
Loaded configurations must have exactly equal up/down counts on every rank.
The model has no pinning field. All physical tensor slices must be populated
and nonzero because the BMPS backend cannot compress zero boundary tensors.
There is no automatic repair of a vanishing total projected amplitude.
Configured `WarmUp` runs even for loaded configurations, since a plain-chain
configuration is not certified thermalized for the projected distribution.

Projected final/lowest outputs use `WavefunctionBase` and carry
`spin_inversion_metadata.txt`; periodic checkpoints inherit the same marker from
`CheckpointBasePath`. Resume with the same parity. Keep this marker with any
copied tensors; when copying a periodic `step_N` directory by itself, also copy
its parent's marker into it. `mc_measure` supports these states when configured with the same
`SpinInversionParity`; plain measurement still rejects marked projected tensors.
Unmarked legacy inputs remain accepted as starting PEPS tensors. A legacy state
which was projected outside this driver must have its parity supplied explicitly.
Both branches currently run sequentially within each MPI chain.
The axis update of section 4.8 cannot be combined with it (`MCAxisUpdate=true` throws).

#### Spatial C4/D4 projection (VMC and measurement)

Both algorithm files accept `PointGroup` (`None`, the default; `C4`; or `D4`),
`PointGroupIrrep`, and the independent `SpinInversionParity` described above.
`C4` contains the four rotations; `D4` adds the four in-plane reflections.
Spin inversion is an independent label operation, not a spatial reflection.

| Group | `PointGroupIrrep` | Meaning |
| --- | --- | --- |
| `None` | `A1` (default) | No spatial projection |
| `C4` | `0` (default), `1`, `2`, `3` as JSON strings | Rotation eigenvalue `i^k` |
| `D4` | `A1` (default), `A2`, `B1`, `B2` | One-dimensional spatial irreps |
| `D4` | `E` | Central projection onto the entire two-dimensional E isotypic subspace |

For D4, `M(r,c)=(r,L-1-c)` defines the reflection: A1/A2 have rotation eigenvalue
+1, B1/B2 have rotation eigenvalue -1, and the suffix 1/2 gives reflection parity
+1/-1. E does not select a rotation eigenvector within its doublet.

Spatial projection requires `Lx=Ly >= 2`, OBC, no removed corners, and
`SquareHeisenberg` or `SquareXY`, and `MCRestrictU1=true` because the projected
exchange updater preserves spin counts. Isotropic J1-J2 couplings are supported, with
J2 on both diagonals. Spin inversion can be independently disabled (`0`) or
assigned parity `+1`/`-1`; when enabled, even site count, `MCRestrictU1=true`,
and Sz=0 loaded configurations on every MPI rank are required. `MCAxisUpdate`
must remain false. Without spatial projection, the existing spin-only VMC
restriction J2=0 remains in effect.

C4 sectors 1 and 3 need complex tensors: configure with `-DRealCode=OFF` and
create/load a complex TPS with that build. Real builds reject these sectors.
D4 and C4 sectors 0/2 also work in the default real build. Tensor file scalar
types must match the executable; changing the build flag does not convert old
wavefunctions.

The coherent amplitude is `sum_g chi(g)* psi(g^-1 x)`, with conjugated
characters and the conventional irrep-dimension/group-order normalization.
All branches share the same base tensors. VMC energy, gradient, SR/MinSR, and
measurement use that amplitude, including interference between branches.
The current reference implementation recomputes canonical contractions for every
proposed configuration; C4/D4 runs are substantially more expensive than the
plain cached sampler. Configured warmup always runs, including for loaded chains.
A starting state annihilated by the chosen projector fails explicitly; choose a
seed with nonzero weight in the desired sector. Increasing contraction accuracy
is still necessary to control BMPS truncation error.

Use matching projection settings for optimization, continuation, and measurement.
Spatial final/lowest and periodic snapshots carry `point_group_metadata.txt`,
including the group, irrep, and spin parity. Keep the marker with copied tensors.
When spin inversion is also enabled, retain its `spin_inversion_metadata.txt`
marker too. A mismatching or plain driver rejects marked tensors; unmarked base
PEPS tensors remain valid starting parameters for projection.

Complete 4x4 D4 A1, spin-even algorithm examples are
`params/quickstart/vmc_cluster_4x4_obc_d4_a1_spin_even_n16.json` and
`params/quickstart/measure_cluster_4x4_obc_d4_a1_spin_even_n16.json`.
Use both with the same 4x4 open-boundary physics file and wavefunction base.

#### 4.8 Axis update (opt-in, OBC)

`MCAxisUpdate=true` adds the PEPS row/column ("axis") autoregressive update to every OBC Monte
Carlo sweep of `vmc_optimize` and `mc_measure`. One sweep then runs
`MCUpdateSquareAxisAutoregressiveOBC`, which redraws every row, then every column, from its
boundary-MPS window and commits the draw without an accept/reject step, followed by the usual
local `MCUpdateSquareTNN3SiteExchangeOBC` sweep (composed with `qlpeps::MCUpdateSequence`).
Consider it when the local updater alone mixes slowly (long autocorrelation times).

Keys (VMC and measurement algorithm files):

- `MCAxisUpdate` (JSON bool, default `false`)
- `MCAxisUpdateDmin` (int, default `1`): `D_min` of the slice MPS
- `MCAxisUpdateDmax` (int, default `Dbmps_max`): `D_max` of the slice MPS (the slice bond
  dimension). Recommended: `Dbmps_max`.
- `MCAxisUpdateTruncErr` (double, default `0.0`): truncation error of the slice MPS

Default off: with `MCAxisUpdate` absent or `false` both programs run exactly as before (the
local updater alone, the same logs). The other three keys are then only type-checked.

Example (entries of the `CaseParams` object of an OBC VMC or measurement algorithm file):

```json
"MCAxisUpdate": true,
"MCAxisUpdateDmax": 16
```

Constraints when `MCAxisUpdate=true` (checked when the parameters are parsed):

- `BoundaryCondition=Open`; `SpinInversionParity=0`.
- `MCAxisUpdateDmax > 0` (its default needs `Dbmps_max`), `1 <= MCAxisUpdateDmin <=
  MCAxisUpdateDmax`, `0 <= MCAxisUpdateTruncErr < 1`.
- When the sampling conserves S_z (a `-DU1SYM` build, or `MCRestrictU1=true`):
  `MCAxisUpdateDmax >= floor(max(Lx, Ly) / 2) + 1`, e.g. at least 7 on 12x12. A bond of a row
  or column can carry that many S_z sectors, and the update keeps at least one singular value
  per sector, so a smaller value can make it throw during sampling (PEPS design
  `docs/dev/design/algorithms/axis-autoregressive-sampling.md`, section 5.4).

Dense versus `-DU1SYM` builds:

| Build | `MCRestrictU1` | Axis moves | Combined chain samples |
|---|---|---|---|
| dense (default) | `true` (default) | keep S_z through the count table N_up `{{1},{0}}` (label 0 = up) | the S_z sector of the initial configuration |
| dense | `false` | change S_z (no table) | all S_z sectors; the local 3-site exchange alone keeps S_z |
| `-DU1SYM` | either | keep S_z, which the tensors carry (no table) | the S_z sector of the tensors |

Keep `MCRestrictU1=true` for a dense state that conserves S_z without carrying it, e.g. one
obtained by simple update from a Neel start: without the table the axis update does not protect
that hidden sector. At a small `MCAxisUpdateDmax` the chain can leave it, and even at a moderate
one it can be confined to part of it, usually with no error (PEPS design, section 5.7).

Accuracy: every draw is committed. With an untruncated slice MPS each move is the exact
heat-bath move in its window (up to the BMPS truncation all OBC updaters share); a truncated
slice MPS adds a bias that is not corrected. Scan `MCAxisUpdateDmax` upward: observables should
not move within error bars. There is no Metropolis-Hastings option.

Cost: the axis update runs in addition to the local sweep and costs more than it.

- Each row or column visit samples its slice. The PEPS design (section 4.13) estimates this
  sampling at about two BMPS compressions of the slice (two labels per site).
- The pass also performs its own BMPS compressions: it grows the boundary MPS above and to the
  left of the visited slice and regrows the one to the right. Its row pass consumes the boundary
  MPS below as it moves down, so the following 3-site sweep must regrow it.
- Measured on 12x12, D=8 Heisenberg states with `Dbmps_max = MCAxisUpdateDmax = 32` (PEPS
  `profiler/README.md`, section "Measured", table "Heisenberg 12x12 D = 8 at N_up = 72", workload
  `heisenberg-sector`; Release build, one thread): a sweep with the axis update took 2.9 times
  as long as a 3-site sweep alone with dense tensors and the count table (18.4 s against
  6.44 s), and 3.8 times with `-DU1SYM` tensors (4.79 s against 1.25 s). It performed 77 BMPS
  compressions instead of 33.
- The energy (and gradient) evaluation of a sample does not change, so the cost per sample grows
  less. With one sweep and one energy evaluation per sample (6.2 s dense, 1.0 s `-DU1SYM` in
  that measurement) it grew about 1.9 times (dense) and 2.6 times (`-DU1SYM`). More sweeps
  between samples bring it closer to the sweep ratio.

Compare autocorrelation times per wall-clock time, not per sweep.

Logs with `MCAxisUpdate=true`:

- At start, rank 0 prints the choice, e.g. for a dense build:

  ```text
  [info] MCAxisUpdate=true: each OBC sweep runs MCUpdateSquareAxisAutoregressiveOBC (rejection-free row and column moves), then MCUpdateSquareTNN3SiteExchangeOBC.
  [info]   slice MPS: D_min=1 D_max=16 trunc_err=0
  [info]   count table: N_up {{1},{0}} (dense tensors, MCRestrictU1=true: S_z is conserved)
  [info]   [MC acceptance] components: 0 = changed rows / Ly, 1 = changed columns / Lx, 2 = local 3-site exchange acceptance; [MC updater] lines carry child0.axis.* counts.
  ```

- `[MC acceptance] chains=<ranks> component=<i> mean=... min=... max=...
  zero_rate_fraction=...`, and the optimizer's `Accept rate = [...]`, now have three
  components instead of one. `mean` averages over chains (MPI ranks) each chain's average over
  its sweeps.
  - Component 0: fraction of rows whose labels changed in a sweep; component 1: the same for
    columns. Every draw is accepted, so these measure how much the slices move, not an
    acceptance rate. Values near 0 mean the slices rarely change: a strongly ordered state or a
    stuck chain.
  - Component 2: acceptance of the local 3-site exchange, the only component without the axis
    update.
- `[MC updater] name=child0.axis.<entry> value=<v> reduction=<sum|max>` lines: VMC prints one
  block per energy evaluation (one per iteration, plus line-search or step-selector
  evaluations; the first block also covers the warm-up), `mc_measure` one block at the end,
  after the data are written. Each block covers the sweeps since the previous one.
  - `row_visits`, `rows_changed`, `col_visits`, `cols_changed` (summed over ranks): the
    pooled change fraction of rows is `rows_changed / row_visits`, likewise for columns.
  - `zero_support_old` (sum): changed visits whose old slice had probability exactly 0 in the
    slice MPS. It should be 0. Otherwise: with dense tensors and `MCRestrictU1=false`, set
    `MCRestrictU1=true` (the state probably conserves S_z without carrying it; see the paragraph
    on such states above); else raise `MCAxisUpdateDmax`. A count of 0 does not rule out the
    confinement described there, which leaves that probability small but not exactly 0.
  - `max_discarded_weight`, `max_slice_bond_dim`, `max_block_count` (maxima over ranks): the
    discarded weight of the slice MPS is 0 when untruncated and is a warning signal, not a bias
    estimate; when it is large, raise `MCAxisUpdateDmax`. `max_block_count` is at most
    `floor(max(Lx, Ly) / 2) + 1` when S_z is conserved, 1 otherwise.

### 5) `measure_algorithm_params.json`

Required baseline keys:

- `MC_total_samples`
- `WarmUp`
- `MCLocalUpdateSweepsBetweenSample`

Backend keys:

- same BMPS (OBC) / TRG or HOTRG (PBC) requirement as VMC, selected by
  `BoundaryCondition` and, for PBC, by `PBCContractor`

Optional keys:

- `ThreadNum` (default `1`)
- IO keys (`WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`) with same defaults as VMC
- `MCRestrictU1` and `InitialConfigStrategy`
- axis-update keys `MCAxisUpdate`, `MCAxisUpdateDmin`, `MCAxisUpdateDmax`,
  `MCAxisUpdateTruncErr` (default off; same semantics and constraints as section 4.8)

Runtime effects:

- loads SITPS from `WavefunctionBase + "final"` when available
- uses warm-start logic from `configuration{rank}` files

### 6) Accepted `MPSCompressScheme` values

- `0` or `"SVD"`
- `1` or `"Variational2Site"`
- `2` or `"Variational1Site"`

### 7) Validation and Hard-Fail Conditions

| Condition | Failure behavior |
|---|---|
| Missing `ModelType` in physics | Throws invalid argument (required key) |
| Invalid `BoundaryCondition` text | Throws invalid argument |
| OBC without `Dbmps_max` | Throws invalid argument (`OBC requested but BMPS params are missing`) |
| PBC without TRG required keys | Throws invalid argument (`TRGDmin`, `TRGDmax`, `TRGTruncErr` required) |
| PBC with `PBCContractor=HOTRG` but without HOTRG required keys | Throws invalid argument (`HOTRGDmin`, `HOTRGDmax`, `HOTRGTruncErr` required) |
| Invalid `PBCContractor` text | Throws invalid argument (`PBCContractor must be TRG or HOTRG.`) |
| Missing `MC_total_samples` / `WarmUp` / `MCLocalUpdateSweepsBetweenSample` | Parse failure in MC param parser |
| SR optimizer without CG required keys | Parse failure in enhanced optimizer parser |
| Step selectors used with non-SGD/SR optimizer | Throws invalid argument |
| Step selectors combined with LR scheduler | Throws invalid argument |
| Selector numeric constraints violated | Throws invalid argument |
| Active advanced-stop + `AdvancedStopEnergyAbsTol` / `AdvancedStopEnergyRelTol` / `AdvancedStopLambdaRelTol` < 0 | Throws invalid argument |
| Active advanced-stop + `AdvancedStopPatience <= 0` or `AdvancedStopMinSteps <= 0` | Throws invalid argument |
| `loop_update` with `ModelType != SquareHeisenberg` | Throws invalid argument |
| `loop_update` with `J2 != 0` | Throws invalid argument |
| `loop_update` with non-positive `ArnoldiTol` / `LoopInvTol` / `FETTolerance` / `CGTol` | Throws invalid argument |
| `loop_update` with non-positive `ArnoldiMaxIter` / `FETMaxIter` / `CGMaxIter` / `CGResidueRestart` | Throws invalid argument |
| `loop_update` with `CGDiagShift < 0` | Throws invalid argument |
| `TauScheduleEnabled=true` but `TauScheduleTaus` missing/empty | Throws invalid argument |
| `TauScheduleTaus` contains non-positive or malformed item | Throws invalid argument |
| `TauScheduleStepCaps` contains non-positive or malformed item | Throws invalid argument |
| `TauScheduleStepCaps` count != `TauScheduleTaus` count | Throws invalid argument |
| `TauScheduleEnabled=true`, `TauScheduleStepCaps` omitted, and `Step <= 0` | Throws invalid argument |
| `TauScheduleRequireConverged=true` but advanced-stop is disabled | Throws invalid argument |
| `TauScheduleDumpDir` empty/whitespace-only | Throws invalid argument |
| `LBFGSStepMode` invalid | Throws invalid argument |
| `LBFGSHistorySize == 0` | Throws invalid argument |
| Strong-Wolfe inequalities violated | Throws invalid argument |
| SITPS boundary != physics boundary (VMC/measure load) | Runtime error and program exit |
| `MCAxisUpdate=true` with `BoundaryCondition=Periodic` or `SpinInversionParity != 0` | Throws invalid argument |
| `MCAxisUpdate=true` with `MCAxisUpdateDmax = 0` (also by default without `Dbmps_max`), `MCAxisUpdateDmin` outside `[1, MCAxisUpdateDmax]`, or `MCAxisUpdateTruncErr` outside `[0, 1)` | Throws invalid argument |
| `MCAxisUpdate=true`, S_z conserved (`-DU1SYM` or `MCRestrictU1=true`), `MCAxisUpdateDmax < floor(max(Lx, Ly) / 2) + 1` | Throws invalid argument (the message states the bound) |
| `MCAxisUpdate*` key of the wrong JSON type (e.g. `"true"` as a string) | Throws invalid argument, also with `MCAxisUpdate=false` |

### 8) High-impact runtime notes

- Programs do not auto-fallback from `tpsfinal/` to `tpslowest/`.
- For `Neel` and `ThreeSublatticePolarizedSeed`, odd `Lx*Ly` in fallback path throws.
- `ThreeSublatticePolarizedSeed` may fall back to `Random` on infeasible small lattices with warning.

For failure diagnosis and exact error text examples, see `tutorials/05-troubleshooting.md`.
