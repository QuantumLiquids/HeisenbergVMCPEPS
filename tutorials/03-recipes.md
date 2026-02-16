## Recipes (Copy, Run, Then Tune)

All commands assume you run inside your build directory (for example `build/` or `build_llvm_sdk/`).

### Recipe A (canonical): Square Heisenberg OBC end-to-end

Use this first if you are new to the project.

#### Physics

```json
{
  "CaseParams": {
    "Lx": 4,
    "Ly": 4,
    "J2": 0.0,
    "ModelType": "SquareHeisenberg",
    "BoundaryCondition": "Open"
  }
}
```

#### Simple Update Algorithm

```json
{
  "CaseParams": {
    "TruncErr": 1e-8,
    "Dmin": 4,
    "Dmax": 4,
    "Tau": 0.2,
    "Step": 100,
    "ThreadNum": 1
  }
}
```

Advanced-stop variant (optional):

```json
{
  "CaseParams": {
    "TruncErr": 1e-8,
    "Dmin": 4,
    "Dmax": 4,
    "Tau": 0.2,
    "Step": 400,
    "ThreadNum": 1,
    "AdvancedStopEnabled": true,
    "AdvancedStopEnergyAbsTol": 1e-8,
    "AdvancedStopEnergyRelTol": 1e-10,
    "AdvancedStopLambdaRelTol": 1e-6,
    "AdvancedStopPatience": 3,
    "AdvancedStopMinSteps": 10
  }
}
```

#### VMC Algorithm (SR + BMPS)

```json
{
  "CaseParams": {
    "ThreadNum": 1,

    "Dbmps_min": 8,
    "Dbmps_max": 8,
    "MPSCompressScheme": "SVD",
    "BMPSTruncErr": 1e-12,

    "MC_total_samples": 3200,
    "WarmUp": 100,
    "MCLocalUpdateSweepsBetweenSample": 1,

    "OptimizerType": "SR",
    "MaxIterations": 30,
    "LearningRate": 0.1,
    "CGMaxIter": 50,
    "CGTol": 1e-6,
    "CGResidueRestart": 10,
    "CGDiagShift": 1e-4,
    "NormalizeUpdate": false
  }
}
```

#### Measure Algorithm (BMPS)

```json
{
  "CaseParams": {
    "ThreadNum": 1,

    "Dbmps_min": 8,
    "Dbmps_max": 8,
    "MPSCompressScheme": "SVD",
    "BMPSTruncErr": 1e-12,

    "MC_total_samples": 6400,
    "WarmUp": 100,
    "MCLocalUpdateSweepsBetweenSample": 1
  }
}
```

#### Run

```bash
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
mpirun -n 1 ./mc_measure ../params/physics_params.json ../params/measure_algorithm_params.json
```

### Recipe B: Square Heisenberg PBC/TRG (when you need periodic boundary)

Use this only after Recipe A is stable.

#### Physics

```json
{
  "CaseParams": {
    "Lx": 4,
    "Ly": 4,
    "J2": 0.3,
    "ModelType": "SquareHeisenberg",
    "BoundaryCondition": "Periodic"
  }
}
```

#### Algorithm requirements

Include TRG keys in VMC and measure algorithm JSON:

- `TRGDmin`
- `TRGDmax`
- `TRGTruncErr`
- optional `TRGInvRelativeEps`

Notes:

- PBC requires SITPS generated consistently for PBC.
- If `tpsfinal/` boundary condition differs from physics JSON, run aborts by design.

### Recipe B0: 2x2 PBC tau schedule (large -> small tau)

Use this when you want staged imaginary-time evolution with decreasing tau:
`0.5 -> 0.2 -> 0.1 -> 0.05 -> 0.02`.

Quickstart files:

- `params/quickstart/physics_local_2x2_pbc.json`
- `params/quickstart/simple_update_local_2x2_pbc_tau_schedule.json`

Run:

```bash
./simple_update \
  ../params/quickstart/physics_local_2x2_pbc.json \
  ../params/quickstart/simple_update_local_2x2_pbc_tau_schedule.json
```

Expected extra schedule logs:

- `=== Tau stage 1/5: tau=0.5, step_cap=8 ===`
- `=== Stage result: converged=..., executed_steps=..., stop_reason=... ===`

Machine-readable outputs:

- `tau_schedule/schedule_summary.json`
- `tau_schedule/schedule_summary.csv`

### Recipe B1: Tile SITPS (PBC/OBC) with `sitps_tile`

Use this after `simple_update` when you want to replicate an existing SITPS into a larger lattice.

PBC example (`2x2 -> 8x8`):

```bash
./sitps_tile \
  --input-dir tpsfinal \
  --output-dir tpsfinal_8x8 \
  --target-ly 8 \
  --target-lx 8 \
  --unit-ly 2 \
  --unit-lx 2
```

Python wrapper (same flags, calls the C++ binary):

```bash
python3 ../scripts/sitps_tile.py \
  --input-dir tpsfinal \
  --output-dir tpsfinal_8x8 \
  --target-ly 8 \
  --target-lx 8 \
  --unit-ly 2 \
  --unit-lx 2
```

Then use matching physics in VMC/measure (`Lx=8, Ly=8, BoundaryCondition=Periodic`) and point `WavefunctionBase` to `tps` with a run directory where `tpsfinal/` is your tiled output.

OBC example (`>=3x3` source only):

```bash
./sitps_tile \
  --input-dir tpsfinal \
  --output-dir tpsfinal_obc_10x12 \
  --target-ly 10 \
  --target-lx 12
```

Important constraints:

- OBC mode rejects `--unit-ly/--unit-lx`.
- PBC mode requires `target` dimensions to be integer multiples of `unit`.
- Tool validates virtual-leg consistency and aborts on mismatch.
- Optional strict cross-check against physics JSON:

```bash
./sitps_tile \
  --input-dir tpsfinal \
  --output-dir tpsfinal_8x8 \
  --target-ly 8 \
  --target-lx 8 \
  --unit-ly 2 \
  --unit-lx 2 \
  --physics-json ../params/physics_8x8_pbc.json
```

### Recipe C: Square XY variant

When to use:

- Same workflow and numerics as square Heisenberg, but XY Hamiltonian.

Change only physics file:

```bash
mpirun -n 1 ./vmc_optimize ../params/physics_square_xy.json ../params/vmc_algorithm_params.json
```

### Recipe D: Triangle Heisenberg variant (OBC only)

When to use:

- Triangle-lattice PEPS path and comparison against DMRG benchmark appendix.

Use triangle physics file:

```bash
./simple_update ../params/physics_triangle_heisenberg.json ../params/simple_update_algorithm_params.json
mpirun -n 1 ./vmc_optimize ../params/physics_triangle_heisenberg.json ../params/vmc_algorithm_params.json
mpirun -n 1 ./mc_measure ../params/physics_triangle_heisenberg.json ../params/measure_algorithm_params.json
```

Note:

- Triangle PBC is not supported.

### Recipe E: Resume from `tpslowest` safely

When to use:

- You want to continue optimization from lowest-energy snapshot instead of latest final snapshot.

```bash
cp -r tpslowest/. tpsfinal/
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
```

This is explicit by design; drivers do not auto-fallback to `tpslowest/`.

### Recipe F: Optimizer variants using provided sample JSONs

Fixed-step L-BFGS:

```bash
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_lbfgs_fixed.json
```

Strong-Wolfe L-BFGS:

```bash
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_lbfgs_strong_wolfe.json
```

SR with step selectors:

```bash
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_sr_step_selector.json
```

### Recipe G: Run + plot trajectory (first-class endpoint)

```bash
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
python3 ../plot/workflow/plot_energy_trajectory.py \
  --csv ./energy/energy_trajectory.csv \
  --out ./energy/energy_trajectory.png \
  --title "VMC trajectory" \
  --ref-energy -9.189 --ref-label "exact"
```

### Recipe H: Machine benchmark case

Use the provided machine test inputs:

```bash
mpirun -n 1 ./vmc_optimize \
  ../machine_test/square12x12D8vmc/physics_params.json \
  ../machine_test/square12x12D8vmc/vmc_algorithm_params.json
```

For DMRG comparison context, see `tutorials/appendix-dmrg-benchmark.md`.
