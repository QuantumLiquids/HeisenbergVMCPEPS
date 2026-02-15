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

Runtime effects:

- `Tau` and `Step` control imaginary-time evolution length/granularity.
- `Dmin`/`Dmax` set SU bond-dimension range.
- Output always targets `tpsfinal/` (SITPS) and `peps/`.

### 3) `vmc_algorithm_params.json`

#### 3.1 Required baseline keys

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

#### 3.2 Backend keys (selected by boundary condition)

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

#### 3.3 IO keys

- `WavefunctionBase` (string, default `"tps"`)
  - load path is `WavefunctionBase + "final"` -> usually `tpsfinal/`
- `ConfigurationLoadDir` (string, default `WavefunctionBase + "final"`)
- `ConfigurationDumpDir` (string, default `WavefunctionBase + "final"`)

Runtime effects:

- `configuration{rank}` is loaded from `ConfigurationLoadDir` when available.
- final configuration is dumped to `ConfigurationDumpDir`.

#### 3.4 Optimizer-specific keys

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

#### 3.5 Step selectors (SGD/SR only)

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

#### 3.6 Misc optional knobs

- `LRScheduler`: `ExponentialDecay`, `CosineAnnealing`, `Plateau`
- `ClipNorm`, `ClipValue`
- spike recovery keys (`Spike*`)
- checkpoint keys (`CheckpointEveryNSteps`, `CheckpointBasePath`)
- `MCRestrictU1` (default `true`, currently informational in unified drivers)
- `InitialConfigStrategy` (default `Random`; accepts `Random`, `Neel`, `ThreeSublatticePolarizedSeed`)

### 4) `measure_algorithm_params.json`

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

### 5) Accepted `MPSCompressScheme` values

- `0` or `"SVD"`
- `1` or `"Variational2Site"`
- `2` or `"Variational1Site"`

### 6) Validation and Hard-Fail Conditions

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
| `LBFGSStepMode` invalid | Throws invalid argument |
| `LBFGSHistorySize == 0` | Throws invalid argument |
| Strong-Wolfe inequalities violated | Throws invalid argument |
| SITPS boundary != physics boundary (VMC/measure load) | Runtime error and program exit |

### 7) High-impact runtime notes

- Programs do not auto-fallback from `tpsfinal/` to `tpslowest/`.
- For `Neel` and `ThreeSublatticePolarizedSeed`, odd `Lx*Ly` in fallback path throws.
- `ThreeSublatticePolarizedSeed` may fall back to `Random` on infeasible small lattices with warning.

For failure diagnosis and exact error text examples, see `tutorials/05-troubleshooting.md`.
