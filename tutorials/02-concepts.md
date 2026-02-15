## Concepts (Core PEPS Workflow)

This repository is a driver layer over `qlpeps`. Keep this model in mind:

- `physics_params.json` tells the code what model to solve.
- `*_algorithm_params.json` tells the code how to solve it (backend, MC, optimizer, IO).

Main workflow:

1. `simple_update` prepares PEPS and SITPS (`tpsfinal/`).
2. `vmc_optimize` loads SITPS and optimizes it.
3. `mc_measure` loads SITPS and measures observables.

### 1) Model / Backend Matrix (strict)

| `ModelType` | `BoundaryCondition` | Contraction backend | Status |
|---|---|---|---|
| `SquareHeisenberg` | `Open` / `OBC` | BMPS | Supported |
| `SquareHeisenberg` | `Periodic` / `PBC` | TRG | Supported |
| `SquareXY` | `Open` / `OBC` | BMPS | Supported |
| `SquareXY` | `Periodic` / `PBC` | TRG | Supported |
| `TriangleHeisenberg` | `Open` / `OBC` | BMPS | Supported |
| `TriangleHeisenberg` | `Periodic` / `PBC` | TRG | Not supported |

### 2) Two-file Parameter Model

Program shape is always:

```bash
./program <physics_params.json> <algorithm_params.json>
```

`physics_params.json`:

- `Lx`, `Ly`, `J2`
- `ModelType` (required)
- `BoundaryCondition` (optional, defaults to open)
- `RemoveCorner` (legacy optional key; ignored by unified square/triangle drivers)

`algorithm_params.json`:

- command-specific numerics
- backend knobs (BMPS or TRG)
- MC/optimizer knobs (where relevant)
- IO overrides (`WavefunctionBase`, `ConfigurationLoadDir`, `ConfigurationDumpDir`)

### 3) Runtime Artifacts and Data Flow

Key directories/files:

- `peps/`: PEPS snapshot written by `simple_update`.
- `tpsfinal/`: canonical SITPS load path for `vmc_optimize` and `mc_measure`.
- `tpslowest/`: best-energy snapshot written by optimizer when available.
- `tpsfinal/configuration{rank}`: per-rank MC configurations (warm start state).
- `energy/energy_trajectory.csv`: human-readable VMC trajectory (canonical plotting input).
- `energy/energy_trajectory`, `energy/energy_err_trajectory`: legacy binary trajectory files.

### 4) State Load / Restart Semantics

`vmc_optimize`:

- Requires `tpsfinal/` by default (`WavefunctionBase + "final"`).
- Does not auto-fallback to `tpslowest/`.
- Enforces boundary-condition consistency between loaded SITPS and physics JSON.

`mc_measure`:

- Prefers SITPS from `tpsfinal/`.
- If SITPS is missing, it may attempt legacy TPS loading and conversion (depends on runtime state and PEPS version).
- Also enforces boundary-condition consistency.

Resume from lowest-energy snapshot:

```bash
cp -r tpslowest/. tpsfinal/
```

### 5) MC Configuration Warm Start

Per-rank config load rule:

1. Try `ConfigurationLoadDir/configuration{rank}`.
2. If loaded, that rank is warmed up.
3. If missing, use `InitialConfigStrategy` fallback and perform warmup sweeps.

Default directories:

- `ConfigurationLoadDir` defaults to `WavefunctionBase + "final"`.
- `ConfigurationDumpDir` defaults to `WavefunctionBase + "final"`.

### 6) Measurement Output Contract (current registry output)

Preferred output layout:

- `<measurement_data_dump_path>/stats/*.csv`
- `<measurement_data_dump_path>/samples/psi.csv`

CSV formats:

1. Flat observables: `stats/<key>.csv` with header `index,mean,stderr`
2. Matrix observables: `stats/<key>_mean.csv` and `stats/<key>_stderr.csv` without header

Current square Heisenberg workflows always include:

- `stats/energy.csv`
- `stats/spin_z_mean.csv`, `stats/spin_z_stderr.csv`
- `stats/bond_energy_h_mean.csv`, `stats/bond_energy_h_stderr.csv`
- `stats/bond_energy_v_mean.csv`, `stats/bond_energy_v_stderr.csv`

Shape rules:

- `spin_z_*`: `Ly x Lx`
- OBC: `bond_energy_h_*` is `Ly x (Lx-1)`, `bond_energy_v_*` is `(Ly-1) x Lx`
- PBC: both bond-energy matrices are `Ly x Lx` (periodic wrap included)

Model/backend-dependent optional examples:

- `stats/SmSp_row.csv`
- `stats/SpSm_row.csv`
- `stats/SzSz_all2all.csv`
- `stats/SzSz_all2all_index_map.txt`


### 7) Deprecation Status

- `src/kagome*` drivers are deprecated and not built by default.
- `plot/` contains mixed legacy/experimental scripts.
- Canonical plotting workflow for tutorials is now:
  - `plot/workflow/plot_energy_trajectory.py`

### 8)  Parameters

See `tutorials/04-parameter-reference.md` for exact defaults and validation rules.
