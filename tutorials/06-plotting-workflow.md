## Plotting Workflow (Python-first)

This tutorial defines one canonical plotting path for the PEPS workflow.

### Goal

Generate a publication-friendly optimization trajectory figure from:

- `energy/energy_trajectory.csv`

with:

- energy curve
- uncertainty band (`energy_error`)
- optional gradient-norm panel

### Input Contract

Required CSV header columns:

- `iteration`
- `energy`
- `energy_error`
- `gradient_norm`

Extra columns are ignored.

`energy_error` may be `nan` when the Monte Carlo standard error (MCSE) is
unavailable. The energy and gradient curves still render, with a gap in the
uncertainty band at those points; a missing MCSE does not mean zero uncertainty.
All other values must be finite, and infinities are rejected in every column.

Canonical source file (after VMC):

- `./energy/energy_trajectory.csv`

### Canonical Command

Run from build directory:

```bash
python3 ../plot/workflow/plot_energy_trajectory.py \
  --csv ./energy/energy_trajectory.csv \
  --out ./energy/energy_trajectory.png \
  --title "Square Heisenberg VMC"
```

Add a horizontal reference line (e.g. exact ground-state energy):

```bash
python3 ../plot/workflow/plot_energy_trajectory.py \
  --csv ./energy/energy_trajectory.csv \
  --out ./energy/energy_trajectory.png \
  --title "Square Heisenberg 4x4 OBC D=4" \
  --ref-energy -9.189 --ref-label "exact"
```

Show interactively in addition to saving:

```bash
python3 ../plot/workflow/plot_energy_trajectory.py \
  --csv ./energy/energy_trajectory.csv \
  --out ./energy/energy_trajectory.png \
  --show
```

### Artifact Check

```bash
ls -l ./energy/energy_trajectory.png
```

### Fallback Note for Legacy Binary Trajectories

Some runs may produce legacy binary files:

- `energy/energy_trajectory`
- `energy/energy_err_trajectory`

The canonical script does not parse those binary files directly. If CSV is missing, regenerate by re-running `vmc_optimize` in a build that writes `energy/energy_trajectory.csv`.

### Legacy MATLAB Scripts

Historical `.m` scripts under `plot/` are still kept for project history and custom workflows.

- They are not the canonical tutorial path.
- They may assume project-specific binary data layouts.

Use `plot/workflow/` for stable tutorial-driven plotting.
