## Concepts (High-Level View)

This repository is a thin, user-facing driver layer on top of `qlpeps`. The two things to keep in your head:
- The physics file says *what* you simulate.
- The algorithm file says *how* you simulate it (including which contraction backend is used).

### What happens in a run (pipeline)

Typical workflow:

1) `simple_update`
   - Input: physics + simple-update algorithm params
   - Output: `peps/` and `tpsfinal/`
2) `vmc_optimize`
   - Input: physics + VMC algorithm params, reads `tpsfinal/`
   - Output: updated `tpsfinal/` (and optional `tpslowest/`, energy CSVs)
3) `mc_measure`
   - Input: physics + measurement algorithm params, reads `tpsfinal/`
   - Output: measurement dumps

### Two-file parameter model (what goes where)

Drivers are called as:
```bash
./prog <physics_params.json> <algorithm_params.json>
```

- `physics_params.json` (model definition)
  - lattice size (`Lx`, `Ly`)
  - couplings (`J2`, ...)
  - model family (`ModelType`)
  - boundary condition (`BoundaryCondition`)
- `*_algorithm_params.json` (numerics + IO)
  - The structure depends on which program you run.

Simple update (`simple_update_algorithm_params.json`):
- Threads: `ThreadNum`
- Simple-update numerics: `Tau`, `Step`, `Dmin`, `Dmax`, `TruncErr`
- No MC / optimizer parameters.

VMC (`vmc_algorithm_params.json`):
- Threads: `ThreadNum`
- Contraction backend params (selected by `BoundaryCondition`):
  - OBC/BMPS: `Dbmps_*`, `MPSCompressScheme`, optional BMPS variational knobs
  - PBC/TRG: `TRG*`
- Monte Carlo params: `MC_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`
- Optimizer params: `OptimizerType`, `MaxIterations`, `LearningRate`, and method-specific keys (SR/Adam/SGD/AdaGrad)
- Optional IO overrides: `WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`

Measurement (`measure_algorithm_params.json`):
- Threads: `ThreadNum`
- Contraction backend params (BMPS or TRG, same rule as VMC)
- Monte Carlo params: `MC_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`
- Optional IO overrides: `WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`
- No optimizer parameters.

This is a taxonomy, not a priority order. Many keys are equally “first-class”; they just control different
subsystems.

### Parameter taxonomy (mental map)

Think in subsystems, not files:

Physics:
- `ModelType` (required), `Lx`, `Ly`, `J2`, `BoundaryCondition`

Wavefunction IO:
- where to load/save tensors: `WavefunctionBase` -> `tpsfinal/`
- where to load/save per-rank MC configuration: `ConfigurationLoadDir`, `ConfigurationDumpDir`

Contraction backend (chosen by physics boundary condition):
- OBC -> BMPS contraction parameters (`Dbmps_*`, `MPSCompressScheme`, optional BMPS variational knobs)
- PBC (square only) -> TRG contraction parameters (`TRG*`)

Monte Carlo:
- `MC_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`

Optimizer (VMC only):
- `OptimizerType`, `MaxIterations`, `LearningRate`, and SR/Adam/SGD specifics

### Boundary condition is a branch, but not the only branch

Boundary condition is one important branch because it selects a different contraction backend:

```json
{
  "CaseParams": {
    "BoundaryCondition": "Open"
  }
}
```

Supported values: `Open` / `OBC` / `Periodic` / `PBC` (case-insensitive).

Backend selection and strict consistency:
- OBC (`Open`) -> BMPS backend
  - requires BMPS parameters in algorithm JSON
- PBC (`Periodic`, square only) -> TRG backend
  - requires TRG parameters in algorithm JSON
- The loaded `tpsfinal/` must have the same boundary condition as `physics_params.json` (enforced).

### Contraction parameters (required/optional)

OBC (BMPS):
- Required: `Dbmps_max`
- Optional: `Dbmps_min` (default = `Dbmps_max`)
- Optional: `MPSCompressScheme` (default = `SVD`)
- Optional: `BMPSTruncErr` (default = 0.0; fallback key `TruncErr`)
- Optional (variational compression only): `BMPSConvergenceTol` (default = 1e-12 if not set), `BMPSIterMax` (default = 10)

PBC (TRG):
- Required: `TRGDmin`, `TRGDmax`, `TRGTruncErr`
- Optional: `TRGInvRelativeEps` (default = 1e-12)

### MPSCompressScheme is string-friendly

You can write:
- `0` or `"SVD"`
- `1` or `"Variational2Site"`
- `2` or `"Variational1Site"`

### Model dispatch (another branch)

`ModelType` selects solver logic:
- `"SquareHeisenberg"`
- `"SquareXY"`
- `"TriangleHeisenberg"` (OBC only)

For square PBC, the solver supports up to J1-J2 (NN + NNN).
