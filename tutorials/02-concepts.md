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
- Monte Carlo params: `MC_total_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`
- Optimizer params: `OptimizerType`, `MaxIterations`, `LearningRate`, and method-specific keys (SR/Adam/SGD/AdaGrad/LBFGS)
- Optional IO overrides: `WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`

Measurement (`measure_algorithm_params.json`):
- Threads: `ThreadNum`
- Contraction backend params (BMPS or TRG, same rule as VMC)
- Monte Carlo params: `MC_total_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`
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
- `MC_total_samples`, `WarmUp`, `MCLocalUpdateSweepsBetweenSample`

Optimizer (VMC only):
- `OptimizerType`, `MaxIterations`, `LearningRate`, and SR/Adam/SGD/AdaGrad/LBFGS specifics
- Optional iterative step selectors are available for SGD/SR (`InitialStepSelector*`, `AutoStepSelector*`) and are not used by LBFGS.

### Monte Carlo configuration (warm start vs cold start)

During VMC optimization and measurement, each MPI rank maintains a Monte Carlo
*configuration* — an assignment of spin-up / spin-down to every lattice site.
The quality of this configuration matters: a good configuration (one that is
already typical for the current wavefunction) lets the sampler produce useful
samples immediately, while a random configuration needs many "warmup" sweeps
before the Markov chain equilibrates.

**How load-from-disk works.** Both `vmc_optimize` and `mc_measure` try to
load a previously saved configuration before they start sampling:

1. Look for the file `ConfigurationLoadDir/configuration{rank}`
   (e.g. `tpsfinal/configuration0` for MPI rank 0).
2. If the file exists and its shape matches the lattice, load it and
   mark this rank as **warmed up** — warmup sweeps are skipped.
3. If the file is missing (e.g. first run, or you added new MPI ranks),
   initialize according to `InitialConfigStrategy` and run the full
   `WarmUp` sweeps specified in the algorithm JSON.

**InitialConfigStrategy (fallback only).**
- `Random` (default): half-up / half-down random occupancy.
- `Neel`: checkerboard AFM pattern with randomly selected phase (two
  sublattice-to-spin mappings are sampled uniformly).
- `ThreeSublatticePolarizedSeed`: 3-sublattice polarized seed.
  - A sublattice: all up.
  - B/C sublattices: choose up spins to be as close as possible to each
    sublattice's 1/4 occupation while enforcing exact total `Sz=0` (for even `Lx*Ly`).
- This strategy is only used when `configuration{rank}` cannot be loaded.
- `Neel` requires even `Lx*Ly` in this fallback path; odd site count throws
  an explicit runtime error.
- For `Neel` with periodic boundary conditions, a true unfrustrated
  checkerboard AFM requires both `Lx` and `Ly` even. If either is odd,
  wrap bonds are frustrated; code prints a per-rank warning and still uses
  the checkerboard seed.
- `ThreeSublatticePolarizedSeed` also requires even `Lx*Ly` in this fallback
  path; odd site count throws an explicit runtime error.
- If `ThreeSublatticePolarizedSeed` is infeasible on a small-size lattice,
  code falls back to `Random` and prints per-rank warning.
- `ThreeSublatticePolarizedSeed` can be used with non-triangle models, but a
  per-rank warning is printed.

After a run finishes, each rank saves its final configuration to
`ConfigurationDumpDir/configuration{rank}`, so the next run can pick up
where this one left off.

**Default directories.** If you do not set `ConfigurationLoadDir` /
`ConfigurationDumpDir` in the algorithm JSON, both default to
`WavefunctionBase + "final"` — which is normally `tpsfinal/`. This means
configurations live alongside the wavefunction tensors, and a simple
re-run of `vmc_optimize` or `mc_measure` automatically reuses them.

**Scaling to more ranks.** If you restart a job with more MPI ranks than
before, the existing ranks load their saved configurations (warm start),
while the new ranks start from random (cold start with warmup). This is
handled automatically — no manual intervention needed.

**Resuming from the lowest-energy snapshot.** The optimizer may save a
`tpslowest/` snapshot. To resume from it, copy `tpslowest/*` into
`tpsfinal/` (including any `configuration{rank}` files that are there).
The programs only look in `tpsfinal/` by default; they never auto-fallback
to `tpslowest/`.

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
