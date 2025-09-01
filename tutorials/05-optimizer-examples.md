## Optimizer parameter examples

These snippets plug directly into your VMC algorithm JSON (e.g., `debug/vmc_quick.json` under `CaseParams`). Only set the keys you need – unspecified fields use sensible defaults, parsed by `EnhancedVMCUpdateParams`.

### Common fields
- OptimizerType: "StochasticReconfiguration" (或简写 "SR") | "SGD" | "Adam" | "AdaGrad"
- MaxIterations, LearningRate
- LRScheduler: "ExponentialDecay" | "CosineAnnealing" | "Plateau"（可选，不填=无调度器）
- ClipNorm, ClipValue（可选，不填=不裁剪；仅对一阶优化器生效）
- WavefunctionBase: "tps" → uses `tpsfinal/` and `tpslowest/`
- ConfigurationLoadDir/ConfigurationDumpDir: directory with per-rank `configuration{rank}` (default `tpsfinal/`)

### BMPS 截断参数的默认
- `Dbmps_max` 必填
- `Dbmps_min` 可选，未给出时默认等于 `Dbmps_max`

### SGD 默认
- `Momentum` 可选，未给出视为 0
- `Nesterov` 可选，未给出视为 false

### Stochastic Reconfiguration (SR)
支持两种写法：
```json
{ "OptimizerType": "SR", "MaxIterations": 5, "LearningRate": 0.1 }
```
或完整写法：
```json
{
  "OptimizerType": "StochasticReconfiguration",
  "MaxIterations": 5,
  "LearningRate": 0.1,
  "CGMaxIter": 50,
  "CGTol": 1e-6,
  "CGResidueRestart": 10,
  "CGDiagShift": 0.01,
  "NormalizeUpdate": false
}
```

### SGD (with Momentum / Nesterov)
```json
{
  "OptimizerType": "SGD",
  "MaxIterations": 10,
  "LearningRate": 0.05,
  "Momentum": 0.9,
  "Nesterov": true,
  "WeightDecay": 0.0,
  "ClipNorm": 1.0,
  "ClipValue": 0.1 
}
```

### Adam
```json
{
  "OptimizerType": "Adam",
  "MaxIterations": 10,
  "LearningRate": 0.01,
  "Beta1": 0.9,
  "Beta2": 0.999,
  "Epsilon": 1e-8,
  "WeightDecay": 0.0,
  "LRScheduler": "CosineAnnealing",
  "MinLearningRate": 1e-4
}
```

### AdaGrad
```json
{
  "OptimizerType": "AdaGrad",
  "MaxIterations": 10,
  "LearningRate": 0.05,
  "Epsilon": 1e-8,
  "InitialAccumulator": 0.0
}
```

### Learning-rate schedulers
- ExponentialDecay:
```json
{ "LRScheduler": "ExponentialDecay", "DecayRate": 0.95, "DecaySteps": 10 }
```
- CosineAnnealing:
```json
{ "LRScheduler": "CosineAnnealing", "MinLearningRate": 1e-4 }
```
- Plateau (energy-based):
```json
{ "LRScheduler": "Plateau", "PlateauFactor": 0.5, "PlateauPatience": 10, "PlateauThreshold": 1e-4 }
```

### Minimal end-to-end quick config
```json
{
  "TruncErr": 1e-10,
  "ThreadNum": 1,
  "MC_samples": 10,
  "WarmUp": 5,
  "MCLocalUpdateSweepsBetweenSample": 1,
  "Dbmps_min": 4,
  "Dbmps_max": 4,
  "MPSCompressScheme": 1,
  "OptimizerType": "StochasticReconfiguration",
  "MaxIterations": 1,
  "LearningRate": 0.1,
  "WavefunctionBase": "tps",
  "ConfigurationLoadDir": "",
  "ConfigurationDumpDir": ""
}
```

