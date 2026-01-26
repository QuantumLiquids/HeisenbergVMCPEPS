## Workflows (short)

### Simple Update (prepare state)
```bash
./simple_update params/physics_params.json params/simple_update_algorithm_params.json
```
Outputs: `peps/`, and dumps a split-index TPS to `tpsfinal/` 用于 VMC。

### VMC Optimization (SR)
```bash
mpirun -n 1 ./vmc_optimize params/physics_params.json params/vmc_algorithm_params.json
```
Outputs: `peps/`, `tpsfinal/`（优化后 SITPS），以及 `energy/energy_trajectory.csv`（CSV）。
实现要点：入口使用上层 API `qlpeps::VmcOptimize(...)`，简化模板显式参数。
示例（伪代码）：
```
auto params = CreateVMCOptimizerParams(...);
SplitIndexTPS sitps; sitps.Load("tpsfinal/");
auto exec = qlpeps::VmcOptimize(params, sitps, MPI_COMM_WORLD, SquareSpinOneHalfXXZModel{}, MCUpdateSquareTNN3SiteExchange{});
```

Format for `energy/energy_trajectory.csv`:
```
update_index  energy_per_site  stderr
```
If different on real runs, replace with actual schema here.

### Monte Carlo Measurement
```bash
./mc_measure src/physics_params.json src/measure_algorithm_params.json
```
Outputs (placeholders):
- `data/energy.dat` → columns: sample_index, energy_per_site
- `data/correlation.dat` → schema TBD; will be updated after validation

### Data flow and load order
- Simple Update writes: `tpsfinal/` (SplitIndexTPS) and `peps/`.
- VMC wavefunction basename: `tps`（即 tpsfinal/, tpslowest/）。
- VMC load rule: 仅从 `tpsfinal/` 加载。若不存在则报错退出。
  - 如需从“最低能量快照”恢复，请手动将 `tpslowest/` 的内容拷贝到 `tpsfinal/`。
- VMC save: 由优化器内部写 `tpsfinal/`（`tpslowest/` 将在接入能量追踪后更新）。
- Measurement: 读取与 VMC 相同的 SplitIndexTPS，推荐使用 `tpsfinal/`。

### Configuration warmup（目录型配置，按 rank）
- 用法：
  - 在 VMC 参数中可提供 `ConfigurationLoadDir` 和 `ConfigurationDumpDir`（留空则默认到 `tpsfinal/`）
  - 每个进程从 `ConfigurationLoadDir/configuration{rank}` 加载；成功则 warmed_up=true；失败（例如新扩容的 rank）则随机并 warmed_up=false
  - 优化器在 `ConfigurationDumpDir/` 下按 `configuration{rank}` 写回最终配置

Example run (small, for minutes):
```bash
# 1) Simple Update to create sitps_final/
./simple_update params/physics_params.json params/simple_update_algorithm_params.json

# 2) VMC (SR by default; supports SGD(momentum)/AdaGrad/clip via params)
#    可选参数：
#    - WavefunctionBase: "tps"（输出到 tpsfinal/tpslowest）
#    - ConfigurationLoadDir: ""（留空 => 使用 tpsfinal/ 作为目录，从 tpsfinal/configuration{rank} 读取）
#    - ConfigurationDumpDir: ""（留空 => 写回到 tpsfinal/configuration{rank}）
./vmc_optimize params/physics_params.json debug/vmc_quick.json

# 3) Measurement (optional)
./mc_measure params/physics_params.json params/measure_algorithm_params.json
```
