## Parameter Reference

This is a compact reference. For runnable examples, use `03-recipes.md`.

### Physics file: `physics_params.json`

Required:
- `Lx` (int)
- `Ly` (int)
- `ModelType` (string): `SquareHeisenberg`, `SquareXY`, `TriangleHeisenberg`
- `J2` (double) (square: NNN coupling; use 0.0 for J1-only)
- `RemoveCorner` (bool) (legacy; kept for compatibility)

Optional:
- `BoundaryCondition` (string): `Open`/`OBC` (default) or `Periodic`/`PBC`

### Algorithm files (shared keys)

Common:
- `ThreadNum` (int): tensor manipulation threads

IO (VMC / measurement):
- `WavefunctionBase` (string): default `"tps"` (loads `tpsfinal/`)
- `ConfigurationLoadDir` (string): directory containing per-rank `configuration{rank}`
- `ConfigurationDumpDir` (string): directory to write per-rank `configuration{rank}`

Monte Carlo:
- `MC_total_samples` (int): required
- `WarmUp` (int): required
- `MCLocalUpdateSweepsBetweenSample` (int): required
- `MCRestrictU1` (bool): optional, default true (currently informational in unified drivers)
- `InitialConfigStrategy` (string): optional, default `Random`
  - `Random`: half-up / half-down random initial configuration when load fails
  - `Neel`: checkerboard AFM initial configuration with random phase when load fails
  - Note: `Neel` requires even `Lx*Ly` when no per-rank configuration is loaded

Optimizer (VMC only):
- `OptimizerType` (string): `SR`/`StochasticReconfiguration`, `SGD`, `Adam`, `AdaGrad`, `LBFGS` (`L-BFGS` alias accepted)
- `MaxIterations` (int)
- `LearningRate` (double)

SGD-only:
- `Momentum` (double; default 0.0)
- `Nesterov` (bool; default false)
- `WeightDecay` (double; default 0.0)

Adam-only:
- `Beta1` (double; default 0.9)
- `Beta2` (double; default 0.999)
- `Epsilon` (double; default 1e-8)
- `WeightDecay` (double; default 0.0)

AdaGrad-only:
- `Epsilon` (double; default 1e-8)
- `InitialAccumulator` (double; default 0.0)

SR-only:
- `CGMaxIter` (int)
- `CGTol` (double)
- `CGResidueRestart` (int)
- `CGDiagShift` (double, default 0.0)
- `NormalizeUpdate` (bool, default false)

LBFGS-only:
- `LBFGSHistorySize` (int, default 10)
- `LBFGSToleranceGrad` (double, default 1e-5)
- `LBFGSToleranceChange` (double, default 1e-9)
- `LBFGSMaxEval` (int, default 20)
- `LBFGSStepMode` (string, default `Fixed`; accepted `Fixed`, `StrongWolfe`, `kFixed`, `kStrongWolfe`)
- `LBFGSWolfeC1` (double, default 1e-4)
- `LBFGSWolfeC2` (double, default 0.9)
- `LBFGSMinStep` (double, default 1e-8)
- `LBFGSMaxStep` (double, default 1.0)
- `LBFGSMinCurvature` (double, default 1e-12)
- `LBFGSUseDamping` (bool, default true)
- `LBFGSMaxDirectionNorm` (double, default 1e3)
- `LBFGSAllowFallbackToFixedStep` (bool, default false)
- `LBFGSFallbackFixedStepScale` (double, default 0.2)

Iterative step selectors (SGD/SR only):
- `InitialStepSelectorEnabled` (bool, default false)
- `InitialStepSelectorMaxLineSearchSteps` (int, default 3)
- `InitialStepSelectorEnableInDeterministic` (bool, default false)
- `AutoStepSelectorEnabled` (bool, default false)
- `AutoStepSelectorEveryNSteps` (int, default 10)
- `AutoStepSelectorPhaseSwitchRatio` (double, default 0.3)
- `AutoStepSelectorEnableInDeterministic` (bool, default false)

Selector constraints:
- Step selectors are only valid with `OptimizerType` = `SGD` or `SR`.
- Step selectors cannot be used together with `LRScheduler`.
- `LearningRate` must be positive when any selector is enabled.
- `InitialStepSelectorMaxLineSearchSteps > 0` when `InitialStepSelectorEnabled=true`.
- `AutoStepSelectorEveryNSteps > 0` and `AutoStepSelectorPhaseSwitchRatio in [0,1]` when `AutoStepSelectorEnabled=true`.

StrongWolfe validation:
- `LBFGSHistorySize > 0` always.
- In `StrongWolfe` mode: `0 < LBFGSWolfeC1 < LBFGSWolfeC2 < 1`, `LBFGSMaxEval > 0`, `LBFGSToleranceGrad >= 0`.

### Simple update algorithm file: `simple_update_algorithm_params.json`

Required:
- `Tau` (double)
- `Step` (int)
- `Dmin` (int)
- `Dmax` (int)
- `TruncErr` (double): simple-update truncation error (SU only)
- `ThreadNum` (int)

### Contraction backend parameters (choose by boundary condition)

OBC (BMPS) required keys:
- `Dbmps_max` (int)

OBC (BMPS) optional keys:
- `Dbmps_min` (int; default = `Dbmps_max`)
- `BMPSTruncErr` (double; default = 0.0). Backward compatible fallback key: `TruncErr`.
- `MPSCompressScheme` (int or string; default = `"SVD"`)
- `BMPSConvergenceTol` (double; only used for variational compression; default = 1e-12)
- `BMPSIterMax` (int; only used for variational compression; default = 10)

PBC (TRG) required keys:
- `TRGDmin` (int)
- `TRGDmax` (int)
- `TRGTruncErr` (double)

PBC (TRG) optional keys:
- `TRGInvRelativeEps` (double; default = 1e-12)

### MPSCompressScheme values

Accepted values:
- `0` or `"SVD"`
- `1` or `"Variational2Site"`
- `2` or `"Variational1Site"`
