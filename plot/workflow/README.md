## Workflow Plot Scripts (Curated)

This directory contains the supported plotting scripts used by the main tutorials.

### Canonical script

- `plot_energy_trajectory.py`

### CLI contract

```bash
python3 plot/workflow/plot_energy_trajectory.py \
  --csv <energy_csv> \
  --out <output_png> \
  [--title <figure_title>] \
  [--ref-energy <value>] \
  [--ref-label <text>] \
  [--show]
```

Required CSV columns:

- `iteration`
- `energy`
- `energy_error`
- `gradient_norm`

Notes:

- Extra columns are ignored.
- `energy_error` may be `nan` when the Monte Carlo standard error is unavailable;
  the uncertainty band has a gap at those points. All other values must be finite,
  and infinities are rejected in every column.
- Script exits non-zero with actionable error messages on missing/invalid input.
