## Quick Start (OBC, minutes)

This is the shortest end-to-end run on a laptop. It uses OBC (BMPS backend).

### Build
```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
mkdir -p build
cd build
cmake ..
make -j4
```

### Run
```bash
# 1) Simple update: generate peps/ and tpsfinal/
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json

# 2) VMC optimize: reads tpsfinal/
mpirun -n 1 ./vmc_optimize ../params/physics_params.json ../params/vmc_algorithm_params.json

# 3) Measurement (optional): reads tpsfinal/
mpirun -n 1 ./mc_measure ../params/physics_params.json ../params/measure_algorithm_params.json
```

### Outputs (important directories)
- `peps/`: PEPS snapshot for continuing simple update
- `tpsfinal/`: SplitIndexTPS (what VMC/measurement load)
- `tpslowest/`: best snapshot during optimization (if enabled)
- `tpsfinal/configuration{rank}`: per-rank Monte Carlo configuration

### If you want PBC

Square PBC uses TRG contraction and requires TRG parameters.

1) Set in `physics_params.json`:
```json
{ "CaseParams": { "BoundaryCondition": "Periodic" } }
```

2) Provide TRG keys in algorithm JSON:
- `TRGDmin`, `TRGDmax`, `TRGTruncErr` (optional: `TRGInvRelativeEps`)

See `02-concepts.md` and `03-recipes.md`.

