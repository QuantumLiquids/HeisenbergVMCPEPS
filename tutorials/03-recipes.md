## Recipes (Copy/Paste)

All examples assume you run from `build/` and paths are relative to it.

Conventions used below:
- The physics JSON defines the model + boundary condition.
- The algorithm JSON defines numerics. It is grouped in this order to keep files readable:
  1) threads
  2) contraction backend params (BMPS or TRG)
  3) MC params
  4) optimizer params (VMC only)
  5) optional IO overrides

Note: `ModelType` is required in the physics file (no default).

### Recipe A: OBC (BMPS) end-to-end

Physics (OBC):
```json
{
  "CaseParams": {
    "Lx": 4,
    "Ly": 4,
    "J2": 0.0,
    "RemoveCorner": true,
    "ModelType": "SquareHeisenberg",
    "BoundaryCondition": "Open"
  }
}
```

Simple update algorithm:
```json
{
  "CaseParams": {
    "TruncErr": 1e-12,
    "Dmin": 4,
    "Dmax": 4,
    "Tau": 0.2,
    "Step": 10,
    "ThreadNum": 1
  }
}
```

VMC algorithm (BMPS + SR):
```json
{
  "CaseParams": {
    "BMPSTruncErr": 1e-12,
    "ThreadNum": 1,

    "Dbmps_min": 8,
    "Dbmps_max": 8,
    "MPSCompressScheme": "SVD",

    "MC_samples": 20,
    "WarmUp": 10,
    "MCLocalUpdateSweepsBetweenSample": 1,

    "OptimizerType": "SR",
    "MaxIterations": 2,
    "LearningRate": 0.1,
    "CGMaxIter": 50,
    "CGTol": 1e-6,
    "CGResidueRestart": 10,
    "CGDiagShift": 1e-4,
    "NormalizeUpdate": false
  }
}
```

Measure algorithm (BMPS):
```json
{
  "CaseParams": {
    "BMPSTruncErr": 1e-12,
    "ThreadNum": 1,

    "Dbmps_min": 8,
    "Dbmps_max": 8,
    "MPSCompressScheme": "SVD",

    "MC_samples": 50,
    "WarmUp": 20,
    "MCLocalUpdateSweepsBetweenSample": 1
  }
}
```

Run:
```bash
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
mpirun -n 1 ./mc_measure ../params/physics_params.json ../params/measure_algorithm_params.json
```

### Recipe B: PBC (TRG) end-to-end (square only)

Physics (PBC):
```json
{
  "CaseParams": {
    "Lx": 4,
    "Ly": 4,
    "J2": 0.3,
    "RemoveCorner": true,
    "ModelType": "SquareHeisenberg",
    "BoundaryCondition": "Periodic"
  }
}
```

VMC algorithm (TRG + SR):
```json
{
  "CaseParams": {
    "ThreadNum": 1,

    "TRGDmin": 4,
    "TRGDmax": 16,
    "TRGTruncErr": 0.0,
    "TRGInvRelativeEps": 1e-12,

    "MC_samples": 20,
    "WarmUp": 10,
    "MCLocalUpdateSweepsBetweenSample": 1,

    "OptimizerType": "SR",
    "MaxIterations": 2,
    "LearningRate": 0.1,
    "CGMaxIter": 50,
    "CGTol": 1e-6,
    "CGResidueRestart": 10,
    "CGDiagShift": 1e-4,
    "NormalizeUpdate": false
  }
}
```

Measure algorithm (TRG):
```json
{
  "CaseParams": {
    "ThreadNum": 1,

    "TRGDmin": 4,
    "TRGDmax": 16,
    "TRGTruncErr": 0.0,
    "TRGInvRelativeEps": 1e-12,

    "MC_samples": 50,
    "WarmUp": 20,
    "MCLocalUpdateSweepsBetweenSample": 1
  }
}
```

Notes:
- For PBC, make sure the `tpsfinal/` you load was generated with PBC (otherwise the driver will abort).
- Triangle PBC is not supported.

### Recipe C: BMPS variational compression knobs

When `MPSCompressScheme` is variational, you can control convergence:

```json
{
  "CaseParams": {
    "Dbmps_min": 8,
    "Dbmps_max": 16,
    "MPSCompressScheme": "Variational2Site",
    "BMPSConvergenceTol": 1e-10,
    "BMPSIterMax": 20
  }
}
```

### Benchmark recipe: square12x12D8vmc (machine_test)

Run the provided benchmark case (two-file mode):

```bash
cd build
mpirun -n 1 ./vmc_optimize \
  ../machine_test/square12x12D8vmc/physics_params.json \
  ../machine_test/square12x12D8vmc/vmc_algorithm_params.json
```
