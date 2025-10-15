## Parameter reference (short)

### physics_params.json
- Lx, Ly: lattice size (int)
- J2: NNN coupling (double)
- RemoveCorner: for triangular geometry (bool)

### numerical (in algorithm files)
- TruncErr (double), ThreadNum (int)
- Dmin/Dmax (int) 仅用于 Simple Update（控制 PEPS bond dimension），VMC/Measure 不使用

### vmc_algorithm_params.json
- 必填：Dbmps_max
- 可选（有默认值）：
  - Dbmps_min（默认等于 Dbmps_max）
  - MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample（无默认值，需提供，否则程序会报错）
  - MPSCompressScheme（无默认则需提供）
  - WavefunctionBase（默认 "tps" → 目录 `tpsfinal/`, `tpslowest/`）
  - ConfigurationLoadDir, ConfigurationDumpDir（默认均为 `wavefunction_base+"final"`，即 `tpsfinal/`）
  - CGDiagShift（SR 中默认 0.0）
  - NormalizeUpdate（SR 中默认 false）
- 其他：CGMaxIter, CGTol, CGResidueRestart（SR 必填），OptimizerType, MaxIterations, LearningRate
- CGMaxIter, CGTol, CGResidueRestart, CGDiagShift
- ReplicaTest
- OptimizerType, MaxIterations, LearningRate
 - WavefunctionBase (e.g. "tps" -> `tpsfinal/`, `tpslowest/`)
 - ConfigurationLoadDir, ConfigurationDumpDir (目录，文件名按 `configuration{rank}`)
  - 详见: [Optimizer Examples](05-optimizer-examples.md)

### simple_update_algorithm_params.json
- Tau, Step

### measure_algorithm_params.json
- 必填：
  - Dbmps_max
  - MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample
  - MPSCompressScheme
- 可选（有默认值）：
  - Dbmps_min（默认等于 Dbmps_max）
  - WavefunctionBase（默认 "tps" → 目录 `tpsfinal/`, `tpslowest/`）
  - ConfigurationLoadDir, ConfigurationDumpDir（默认均为 `wavefunction_base+"final"`，即 `tpsfinal/`）

### Presets
```json
// fast test (VMC/Measure: Dbmps_* for BMPS, no Dmin/Dmax here)
{"TruncErr":1e-10}

// accurate
{"TruncErr":1e-15}
```
