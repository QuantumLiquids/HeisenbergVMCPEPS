#!/usr/bin/env python3
"""Plot VMC energy trajectory from CSV.

Expected CSV columns:
- iteration
- energy
- energy_error
- gradient_norm
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path

REQUIRED_COLUMNS = ("iteration", "energy", "energy_error", "gradient_norm")


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
      description="Plot VMC energy trajectory from energy_trajectory.csv")
  parser.add_argument("--csv", required=True, help="Input CSV path")
  parser.add_argument("--out", required=True, help="Output PNG path")
  parser.add_argument("--title", default="", help="Optional figure title")
  parser.add_argument("--ref-energy", type=float, default=None,
                      help="Reference energy to draw as horizontal dashed line")
  parser.add_argument("--ref-label", default="exact",
                      help="Label for the reference energy line (default: 'exact')")
  parser.add_argument("--show", action="store_true", help="Show figure window")
  return parser.parse_args()


def fail(message: str, code: int = 2) -> int:
  print(f"error: {message}", file=sys.stderr)
  return code


def parse_float(value: str, column: str, line_no: int) -> float:
  try:
    parsed = float(value)
  except ValueError as exc:
    raise ValueError(
        f"line {line_no}: column '{column}' is not a float: {value!r}") from exc
  if not math.isfinite(parsed):
    raise ValueError(
        f"line {line_no}: column '{column}' is non-finite: {value!r}")
  return parsed


def load_csv(path: Path) -> tuple[list[float], list[float], list[float], list[float]]:
  if not path.exists():
    raise FileNotFoundError(f"CSV file not found: {path}")
  if not path.is_file():
    raise ValueError(f"CSV path is not a regular file: {path}")

  with path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    if not reader.fieldnames:
      raise ValueError("CSV header is missing")

    missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
    if missing:
      raise ValueError(
          "CSV is missing required columns: " + ", ".join(missing))

    iterations: list[float] = []
    energies: list[float] = []
    errors: list[float] = []
    grad_norms: list[float] = []

    for line_no, row in enumerate(reader, start=2):
      if row is None:
        continue
      iteration = parse_float(row["iteration"], "iteration", line_no)
      energy = parse_float(row["energy"], "energy", line_no)
      energy_error = parse_float(row["energy_error"], "energy_error", line_no)
      grad_norm = parse_float(row["gradient_norm"], "gradient_norm", line_no)

      iterations.append(iteration)
      energies.append(energy)
      # Negative error bars are not meaningful for plotting; clamp to zero.
      errors.append(max(0.0, energy_error))
      grad_norms.append(max(0.0, grad_norm))

  if not iterations:
    raise ValueError("CSV has header but no data rows")

  return iterations, energies, errors, grad_norms


def main() -> int:
  args = parse_args()

  csv_path = Path(args.csv)
  out_path = Path(args.out)

  try:
    iterations, energies, errors, grad_norms = load_csv(csv_path)
  except FileNotFoundError as exc:
    return fail(str(exc), code=2)
  except ValueError as exc:
    return fail(str(exc), code=2)

  # Keep matplotlib caches in writable temp location under restricted environments.
  if "MPLCONFIGDIR" not in os.environ:
    mpl_dir = Path("/tmp/matplotlib")
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_dir)

  import matplotlib

  if not args.show:
    matplotlib.use("Agg")

  import matplotlib.pyplot as plt

  fig, axes = plt.subplots(
      nrows=2,
      ncols=1,
      figsize=(9, 7),
      sharex=True,
      gridspec_kw={"height_ratios": [3, 2]})
  ax_energy, ax_grad = axes

  ax_energy.plot(iterations, energies, color="tab:blue", lw=1.8, label="energy")
  if any(err > 0.0 for err in errors):
    lower = [e - err for e, err in zip(energies, errors)]
    upper = [e + err for e, err in zip(energies, errors)]
    ax_energy.fill_between(
        iterations,
        lower,
        upper,
        color="tab:blue",
        alpha=0.2,
        linewidth=0.0,
        label="energy ± error")
  if args.ref_energy is not None:
    ax_energy.axhline(args.ref_energy, color="tab:red", ls="--", lw=1.2,
                      label=f"{args.ref_label} ({args.ref_energy})")
  ax_energy.set_ylabel("Energy")
  ax_energy.grid(True, alpha=0.3)
  ax_energy.legend(loc="best")

  if any(val > 0.0 for val in grad_norms):
    positive_x = [x for x, y in zip(iterations, grad_norms) if y > 0.0]
    positive_y = [y for y in grad_norms if y > 0.0]
    ax_grad.semilogy(positive_x, positive_y, color="tab:orange", lw=1.8,
                     label="gradient norm")
    ax_grad.set_ylabel("Grad Norm (log)")
  else:
    ax_grad.plot(iterations, grad_norms, color="tab:orange", lw=1.8,
                 label="gradient norm")
    ax_grad.set_ylabel("Grad Norm")

  ax_grad.set_xlabel("Iteration")
  ax_grad.grid(True, alpha=0.3)
  ax_grad.legend(loc="best")

  title = args.title.strip() if args.title else "VMC Energy Trajectory"
  fig.suptitle(title)
  fig.tight_layout()

  out_path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(out_path, dpi=160)

  print(f"[ok] plotted {len(iterations)} points to {out_path}")

  if args.show:
    plt.show()
  else:
    plt.close(fig)

  return 0


if __name__ == "__main__":
  sys.exit(main())
