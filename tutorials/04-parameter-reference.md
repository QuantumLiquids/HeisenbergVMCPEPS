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

- Driver runs loop-update sweeps with `LoopUpdatePara` built from this JSON.
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

PBC/TRG:

- required: `TRGDmin`, `TRGDmax`, `TRGTruncErr`
- optional: `TRGInvRelativeEps` (default `1e-12`)

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
- `MCRestrictU1` (default `true`, currently informational in unified drivers)
- `InitialConfigStrategy` (default `Random`; accepts `Random`, `Neel`, `ThreeSublatticePolarizedSeed`)

### 5) `measure_algorithm_params.json`

Required baseline keys:

- `MC_total_samples`
- `WarmUp`
- `MCLocalUpdateSweepsBetweenSample`

Backend keys:

- same BMPS/TRG requirement as VMC, selected by `BoundaryCondition`

Optional keys:

- `ThreadNum` (default `1`)
- IO keys (`WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`) with same defaults as VMC
- `MCRestrictU1` and `InitialConfigStrategy`

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

### 8) High-impact runtime notes

- Programs do not auto-fallback from `tpsfinal/` to `tpslowest/`.
- For `Neel` and `ThreeSublatticePolarizedSeed`, odd `Lx*Ly` in fallback path throws.
- `ThreeSublatticePolarizedSeed` may fall back to `Random` on infeasible small lattices with warning.

For failure diagnosis and exact error text examples, see `tutorials/05-troubleshooting.md`.
