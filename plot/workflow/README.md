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
- Script exits non-zero with actionable error messages on missing/invalid input.
