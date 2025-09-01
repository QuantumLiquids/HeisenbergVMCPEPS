### Machine test: square12x12D8vmc（VMC benchmark）

**模型与编译要求**
- **模型**: square-lattice spin-1/2 Heisenberg（J2=0）
- **对称性**: 必须使用不带 U1 的编译（默认即可）。不要设置 `-DU1SYM=1`。

**目录与数据**
- 目录: `machine_test/square12x12D8vmc/`
- 已提供：
  - `tps/` 与（可能的）`peps/` 初始波函数数据（TPS 为 legacy 版本）
  - `vmc_params.json`, `simple_update_params.json`
  - 运行日志 `.log`


**运行 VMC**
```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/build
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json
# 或使用机器算例自带参数：
mpirun -n 1 ./vmc_optimize ../machine_test/square12x12D8vmc/vmc_params.json ../params/vmc_algorithm_params.json
```

注意：
- VMC 仅从 `tpsfinal/` 加载 SplitIndexTPS；若需从最低能量快照恢复，请手动将 `tpslowest/` 拷贝到 `tpsfinal/`。
- 已提供日志用于速度与收敛的参考比对。


