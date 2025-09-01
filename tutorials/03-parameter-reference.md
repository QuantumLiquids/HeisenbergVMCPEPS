## Parameter reference (short)

### physics_params.json
- Lx, Ly: lattice size (int)
- J2: NNN coupling (double)
- RemoveCorner: for triangular geometry (bool)

### numerical (in algorithm files)
- TruncErr (double), ThreadNum (int)
- Dmin/Dmax (int) 仅用于 Simple Update（控制 PEPS bond dimension），VMC/Measure 不使用

### vmc_algorithm_params.json
- MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample
- Dbmps_min/Dbmps_max, MPSCompressScheme
- CGMaxIter, CGTol, CGResidueRestart, CGDiagShift
- ReplicaTest
- OptimizerType, MaxIterations, LearningRate
 - WavefunctionBase (e.g. "tps" -> `tpsfinal/`, `tpslowest/`)
 - ConfigurationLoadDir, ConfigurationDumpDir (目录，文件名按 `configuration{rank}`)
  - 详见: [Optimizer Examples](05-optimizer-examples.md)

### simple_update_algorithm_params.json
- Tau, Step

### measure_algorithm_params.json
- MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample
- Dbmps_min/Dbmps_max, MPSCompressScheme

### Presets
```json
// fast test (VMC/Measure: Dbmps_* for BMPS, no Dmin/Dmax here)
{"TruncErr":1e-10}

// accurate
{"TruncErr":1e-15}
```
