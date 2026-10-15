## Quickstart Parameter Presets

`local22_pbc` (tiny local run):
- `params/quickstart/physics_local_2x2_pbc.json`
- `params/quickstart/simple_update_local_2x2_pbc.json`
- `params/quickstart/simple_update_local_2x2_pbc_advanced_stop.json` (advanced-stop variant)
- `params/quickstart/simple_update_local_2x2_pbc_tau_schedule.json` (tau-schedule variant)
- `params/quickstart/loop_update_local_2x2_pbc.json`
- `params/quickstart/loop_update_local_2x2_pbc_advanced_stop.json` (advanced-stop variant)
- `params/quickstart/vmc_local_2x2_pbc_robust_n1.json` (default tutorial SR preset)
- `params/quickstart/vmc_local_2x2_pbc_n1.json` (original preset, retained for comparisons)
- `params/quickstart/measure_local_2x2_pbc_n1.json`
- `params/quickstart/vmc_local_2x2_pbc_hotrg_n1.json` (HOTRG contractor variant)
- `params/quickstart/measure_local_2x2_pbc_hotrg_n1.json` (HOTRG contractor variant)

`cluster44_obc` (cluster-oriented run):
- `params/quickstart/physics_cluster_4x4_obc.json`
- `params/quickstart/simple_update_cluster_4x4_obc.json`
- `params/quickstart/vmc_cluster_4x4_obc_n16.json`
- `params/quickstart/measure_cluster_4x4_obc_n16.json`

The workflow script `scripts/run_quickstart_workflow.sh` snapshots the exact
JSON files used into `run/<timestamp>_<profile>/params_used/`, so parameter
changes are fully traceable for later debugging.

### Repeated 2x2 PBC convergence runs

The tutorial and `local22_pbc` workflow use the robust preset with D=2,
eight-sweep simple-update initialization and 1,000 samples per VMC iteration:

```bash
scripts/run_quickstart_workflow.sh \
  --profile local22_pbc \
  --build-dir /absolute/path/to/build
```

The workflow creates a fresh run directory. This matters because `simple_update`
automatically loads an existing `peps/`, and VMC reads and overwrites `tpsfinal/`.
The plotter leaves gaps for unavailable (`nan`) Monte Carlo error bars while
preserving every energy and gradient point. Use `--no-plot` if Matplotlib is not
installed. To compare with the original preset, pass
`--vmc params/quickstart/vmc_local_2x2_pbc_n1.json` explicitly.

Compared with the original preset, the robust preset changes only:

| Parameter | Original | Robust |
|---|---:|---:|
| `CGRelativeTolerance` | `1e-3` | `1e-8` |
| `SRDiagShift` | `1e-4` | `1e-2` |
| `LearningRate` | `0.05` | `0.2` |

The tighter CG tolerance improves the accuracy of the SR linear solve; stronger SR
regularization controls nearly singular tangent directions and permits the
larger step. The historical `CGTol=1e-6` means squared relative residual and is
equivalent to the original `CGRelativeTolerance=1e-3`, not the robust setting.

Validated on 2026-10-15 with HeisenbergVMCPEPS `fb930369`, PEPS-Pro `22d460fd`,
UltraDMRG `7ed427af`, TensorToolkit-Pro `2dc6fe9d`, GCC 14.2, OpenBLAS, real dense
tensors, one MPI process and one BLAS thread. Ten independent 30-step runs all
finished with logged gradient norms below `4e-4` and exact Sz=0 energy errors
below `3e-8` relative to total energy -4. The final ten logged gradient values
decreased in every run. Early MC fluctuations and temporary gradient rises
remain; this is not a guarantee of stepwise monotonicity or bitwise agreement.
The application currently does not expose a fixed Monte Carlo seed in JSON.
Exact energy here refers to the fixed Sz=0 sector sampled by this tutorial.
