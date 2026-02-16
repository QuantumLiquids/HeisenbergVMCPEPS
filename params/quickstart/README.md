## Quickstart Parameter Presets

`local22_pbc` (stable tiny local run):
- `params/quickstart/physics_local_2x2_pbc.json`
- `params/quickstart/simple_update_local_2x2_pbc.json`
- `params/quickstart/simple_update_local_2x2_pbc_advanced_stop.json` (advanced-stop variant)
- `params/quickstart/simple_update_local_2x2_pbc_tau_schedule.json` (tau-schedule variant)
- `params/quickstart/vmc_local_2x2_pbc_n1.json`
- `params/quickstart/measure_local_2x2_pbc_n1.json`

`cluster44_obc` (cluster-oriented run):
- `params/quickstart/physics_cluster_4x4_obc.json`
- `params/quickstart/simple_update_cluster_4x4_obc.json`
- `params/quickstart/vmc_cluster_4x4_obc_n16.json`
- `params/quickstart/measure_cluster_4x4_obc_n16.json`

The workflow script `scripts/run_quickstart_workflow.sh` snapshots the exact
JSON files used into `run/<timestamp>_<profile>/params_used/`, so parameter
changes are fully traceable for later debugging.
