## Two-file parameter system

Keep physics separate from algorithms.

### Files
```
params/
├── physics_params.json
├── simple_update_algorithm_params.json
├── vmc_algorithm_params.json
└── measure_algorithm_params.json
```

### physics_params.json
```json
{
  "CaseParams": { "Lx": 4, "Ly": 4, "J2": 0.0, "RemoveCorner": true }
}
```

### Numerical blocks by algorithm
- Simple Update:
```json
{ "CaseParams": { "TruncErr": 1e-15, "Dmin": 4, "Dmax": 4, "ThreadNum": 2 } }
```
- VMC / Measurement (boundary MPS truncation only):
```json
{ "CaseParams": { "TruncErr": 1e-15, "ThreadNum": 2 } }
```

### VMC add-ons
```json
{
  "CaseParams": {
    "MC_samples": 1000, "WarmUp": 100,
    "MCLocalUpdateSweepsBetweenSample": 1,
    "Dbmps_min": 4, "Dbmps_max": 4, "MPSCompressScheme": 1,
    "CGMaxIter": 100, "CGTol": 1e-8, "CGResidueRestart": 20, "CGDiagShift": 0.01,
    "ReplicaTest": true,
    "OptimizerType": "SR", "MaxIterations": 2, "LearningRate": 0.1
  }
}
```

### Simple Update add-ons
```json
{ "CaseParams": { "Tau": 0.2, "Step": 10 } }
```

### Measurement add-ons
```json
{
  "CaseParams": {
    "MC_samples": 2000, "WarmUp": 1000,
    "MCLocalUpdateSweepsBetweenSample": 2,
    "Dbmps_min": 4, "Dbmps_max": 4, "MPSCompressScheme": 1
  }
}
```

### Usage
```bash
./prog <physics_params.json> <algorithm_params.json>
# example
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
```
