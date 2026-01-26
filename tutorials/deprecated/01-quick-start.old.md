## Quick Start

Follow this to run a tiny end-to-end workflow (minutes on a laptop).

### Build
```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
mkdir -p build && cd build
cmake ..
make -j4
```

### Minimal params (already provided)
- `params/physics_params.json` (4x4, J2=0)
- `params/simple_update_algorithm_params.json` (D=4, Step=10)
- `params/vmc_algorithm_params.json`（或 `debug/vmc_quick.json`）

### Run
```bash
# 1) Prepare state
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json

# 2) Optimize by VMC
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json

# 3) Measure (optional quick check)
./mc_measure ../src/physics_params.json ../src/measure_algorithm_params.json
```

### Outputs（实测）
- `peps/`: 保存 PEPS（用于Simple update续算）
- `tpsfinal/`: 保存 SplitIndexTPS（VMC 唯一加载路径）
- `tpslowest/`: VMC 运行过程中记录的最低能量快照
  - 如需从最低能量快照恢复，请手动将 `tpslowest/` 拷贝到 `tpsfinal/`
- `tpsfinal/configuration{rank}`: 每个进程的 MC 配置
- `energy/energy_trajectory.csv`: 每次迭代的能量与误差（CSV，二进制轨迹已弃用）
