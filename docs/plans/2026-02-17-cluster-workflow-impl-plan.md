# Cluster Workflow Skill + CLI Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the unused MCP server and old skills with a lightweight Python CLI (`scripts/cluster_run.py`) and a Claude Code skill that together automate HPC job setup, submission, monitoring, and result fetching for VMC-PEPS calculations.

**Architecture:** A skill file encodes workflow knowledge and conventions. A Python CLI generates case directories with parameter files and slurm scripts from templates, then submits/monitors/fetches via SSH. No MCP server, no external dependencies beyond Python stdlib.

**Tech Stack:** Python 3.8+ (runtime: stdlib only; dev: pytest), Bash/SSH, Claude Code skills (markdown)

**Design doc:** `docs/plans/2026-02-17-cluster-workflow-skill-cli-design.md`

---

### Task 1: Delete unused infrastructure

Remove the untracked MCP server, old design doc, and old skill that were never deployed.

**Files:**
- Delete: `mcp/run-manager/` (entire directory)
- Delete: `docs/plans/2026-02-12-mcp-run-manager-v01-design.md`
- Delete: `skills/SKILL.md`
- Keep: `skills/susphy_cluster_run_notes.md` (reference for Task 7)

**Step 1: Remove files**

```bash
rm -rf mcp/run-manager/
rm docs/plans/2026-02-12-mcp-run-manager-v01-design.md
rm skills/SKILL.md
```

Verify `mcp/` directory is now empty (or only contains `.gitkeep` if tracked).
Verify `skills/susphy_cluster_run_notes.md` still exists.

**Step 2: Commit**

```bash
git add -A mcp/ docs/plans/2026-02-12-mcp-run-manager-v01-design.md skills/SKILL.md
git commit -m "chore: remove unused MCP run-manager and old skill/design docs"
```

---

### Task 2: CLI scaffold with `new-case` parameter computation

Build the core parameter computation logic that takes CLI inputs and produces
the full set of JSON-ready parameter dictionaries. This is the pure-logic core
with no I/O.

**Files:**
- Create: `scripts/cluster_run.py`
- Create: `tests/test_cluster_run.py`

**Step 1: Write failing tests for parameter computation**

Create `tests/test_cluster_run.py`:

```python
"""Tests for cluster_run.py parameter computation."""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from cluster_run import (
    compute_mc_samples, compute_case_name, build_params,
    parse_sbatch_job_id, parse_lattice_spec,
)


class TestParseSbatchJobId:
    def test_normal_output(self):
        assert parse_sbatch_job_id("Submitted batch job 517588\n") == "517588"

    def test_with_banner_text(self):
        out = "Some MOTD banner\nSubmitted batch job 123456\n"
        assert parse_sbatch_job_id(out) == "123456"

    def test_missing_job_id_raises(self):
        with pytest.raises(ValueError, match="Could not parse"):
            parse_sbatch_job_id("sbatch: error: invalid partition\n")

    def test_empty_output_raises(self):
        with pytest.raises(ValueError, match="Could not parse"):
            parse_sbatch_job_id("")


class TestParseLatticeSpec:
    def test_normal(self):
        assert parse_lattice_spec("8x8") == (8, 8)

    def test_rectangular(self):
        assert parse_lattice_spec("12x8") == (12, 8)

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Invalid lattice"):
            parse_lattice_spec("8by8")

    def test_zero_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            parse_lattice_spec("0x8")


class TestComputeMcSamples:
    def test_small_lattice_clamps_to_minimum(self):
        # 4x4 = 3200 raw, clamped to 10000, rounded to nearest 56
        result = compute_mc_samples(4, 4, 56)
        assert result >= 10000
        assert result % 56 == 0

    def test_medium_lattice_scales(self):
        # 8x8 = 12800 raw, within range, rounded to nearest 56
        result = compute_mc_samples(8, 8, 56)
        assert 10000 <= result <= 20000
        assert result % 56 == 0

    def test_large_lattice_clamps_to_maximum(self):
        # 16x16 = 51200 raw, clamped to 20000, rounded to nearest 56
        result = compute_mc_samples(16, 16, 56)
        assert result <= 20000
        assert result % 56 == 0


class TestComputeCaseName:
    def test_integer_j2(self):
        assert compute_case_name(8, 8, 0.0, 8) == "8x8J2=0D8"

    def test_fractional_j2(self):
        assert compute_case_name(8, 8, 0.5, 8) == "8x8J2=0.5D8"

    def test_different_lattice(self):
        assert compute_case_name(12, 12, 0.0, 5) == "12x12J2=0D5"


class TestBuildParams:
    def test_physics_unit_params(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        pu = p["physics_unit"]
        assert pu["CaseParams"]["Lx"] == 2
        assert pu["CaseParams"]["Ly"] == 2
        assert pu["CaseParams"]["BoundaryCondition"] == "Periodic"

    def test_physics_full_params(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        pf = p["physics_full"]
        assert pf["CaseParams"]["Lx"] == 8
        assert pf["CaseParams"]["Ly"] == 8

    def test_simple_update_has_tau_schedule(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        su = p["simple_update"]["CaseParams"]
        assert su["TauScheduleEnabled"] is True
        assert su["TauScheduleTaus"] == "0.5,0.2,0.1,0.05,0.02"
        assert su["Dmin"] == 8
        assert su["Dmax"] == 8

    def test_vmc_trg_defaults(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        vmc = p["vmc"]["CaseParams"]
        assert vmc["TRGDmin"] == 8
        assert vmc["TRGDmax"] == 16
        assert vmc["OptimizerType"] == "SR"

    def test_measure_doubles_mc_samples(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        vmc_samples = p["vmc"]["CaseParams"]["MC_total_samples"]
        meas_samples = p["measure"]["CaseParams"]["MC_total_samples"]
        assert meas_samples == 2 * vmc_samples

    def test_override_applies(self):
        p = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={"CGMaxIter": "500"},
        )
        assert p["vmc"]["CaseParams"]["CGMaxIter"] == 500

    def test_obc_uses_bmps_not_trg(self):
        p = build_params(
            lx=4, ly=4, d=4, j2=0.0, bc="obc",
            unit_lx=4, unit_ly=4,
            ntasks=16, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        vmc = p["vmc"]["CaseParams"]
        assert "Dbmps_min" in vmc
        assert "TRGDmin" not in vmc
```

**Step 2: Run tests to verify they fail**

```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: `ModuleNotFoundError: No module named 'cluster_run'`

**Step 3: Write the parameter computation functions**

Create `scripts/cluster_run.py` with these functions (no CLI entry point yet):

```python
#!/usr/bin/env python3
"""CLI for generating, submitting, and managing VMC-PEPS cluster runs."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Parameter computation (pure functions, no I/O)
# ---------------------------------------------------------------------------

def parse_sbatch_job_id(stdout: str) -> str:
    """Extract job ID from sbatch output using regex. Raises ValueError on failure."""
    match = re.search(r"Submitted batch job (\d+)", stdout)
    if not match:
        raise ValueError(f"Could not parse sbatch job ID from output: {stdout!r}")
    return match.group(1)


def parse_lattice_spec(spec: str) -> tuple[int, int]:
    """Parse 'LxXLy' string (e.g. '8x8') into (lx, ly). Raises ValueError on bad input."""
    match = re.fullmatch(r"(\d+)[xX](\d+)", spec.strip())
    if not match:
        raise ValueError(f"Invalid lattice spec: {spec!r}. Expected format: 8x8")
    lx, ly = int(match.group(1)), int(match.group(2))
    if lx <= 0 or ly <= 0:
        raise ValueError(f"Lattice dimensions must be positive, got {lx}x{ly}")
    return lx, ly


def compute_mc_samples(lx: int, ly: int, ntasks: int) -> int:
    """Compute MC_total_samples scaled by lattice size, clamped to [10000, 20000]."""
    raw = lx * ly * 200
    clamped = max(10000, min(20000, raw))
    # Round to nearest multiple of ntasks
    return int(math.ceil(clamped / ntasks) * ntasks)


def compute_case_name(lx: int, ly: int, j2: float, d: int) -> str:
    """Build case directory name like 8x8J2=0D8."""
    j2_str = str(j2) if j2 != int(j2) else str(int(j2))
    return f"{lx}x{ly}J2={j2_str}D{d}"


def build_params(
    lx: int,
    ly: int,
    d: int,
    j2: float,
    bc: str,
    unit_lx: int,
    unit_ly: int,
    ntasks: int,
    optimizer: str,
    lr: float,
    iterations: int,
    taus: str,
    step_cap: int,
    overrides: Dict[str, str],
) -> Dict[str, Any]:
    """Build all five parameter dictionaries from CLI inputs."""

    bc_value = "Periodic" if bc == "pbc" else "Open"
    model = "SquareHeisenberg"
    vmc_samples = compute_mc_samples(lx, ly, ntasks)
    meas_samples = 2 * vmc_samples

    physics_unit = {
        "CaseParams": {
            "Lx": unit_lx,
            "Ly": unit_ly,
            "J2": j2,
            "RemoveCorner": True,
            "ModelType": model,
            "BoundaryCondition": bc_value,
        }
    }

    physics_full = {
        "CaseParams": {
            "Lx": lx,
            "Ly": ly,
            "J2": j2,
            "RemoveCorner": True,
            "ModelType": model,
            "BoundaryCondition": bc_value,
        }
    }

    simple_update = {
        "CaseParams": {
            "TruncErr": 1e-8,
            "Dmin": d,
            "Dmax": d,
            "ThreadNum": 1,
            "Tau": float(taus.split(",")[0]),
            "Step": step_cap,
            "TauScheduleEnabled": True,
            "TauScheduleTaus": taus,
            "TauScheduleRequireConverged": True,
            "TauScheduleDumpEachStage": True,
            "TauScheduleDumpDir": "tau_schedule",
            "TauScheduleAbortOnStageFailure": True,
            "AdvancedStopEnabled": True,
            "AdvancedStopEnergyAbsTol": 1e-6,
            "AdvancedStopEnergyRelTol": 1e-8,
            "AdvancedStopLambdaRelTol": 1e-5,
            "AdvancedStopPatience": 3,
            "AdvancedStopMinSteps": 10,
        }
    }

    # VMC params differ for PBC (TRG) vs OBC (BMPS)
    vmc_case: Dict[str, Any] = {
        "ThreadNum": 1,
        "InitialConfigStrategy": "Neel",
        "MC_total_samples": vmc_samples,
        "WarmUp": 100,
        "MCLocalUpdateSweepsBetweenSample": 1,
        "OptimizerType": optimizer,
        "MaxIterations": iterations,
        "LearningRate": lr,
        "CGMaxIter": 200,
        "CGTol": 1e-6,
        "CGResidueRestart": 10,
        "CGDiagShift": 1e-4,
        "NormalizeUpdate": False,
        "InitialStepSelectorEnabled": True,
        "InitialStepSelectorMaxLineSearchSteps": 3,
        "InitialStepSelectorEnableInDeterministic": False,
        "AutoStepSelectorEnabled": True,
        "AutoStepSelectorEveryNSteps": 10,
        "AutoStepSelectorPhaseSwitchRatio": 0.3,
        "AutoStepSelectorEnableInDeterministic": False,
        "CheckpointEveryNSteps": 10,
        "CheckpointBasePath": "checkpoints",
        "SpikeAutoRecover": True,
        "SpikeMaxRetries": 2,
        "SpikeFactorErr": 100.0,
        "SpikeFactorGrad": 1e10,
        "SpikeFactorNGrad": 10.0,
        "SpikeSRMinItersSuspicious": 1,
        "SpikeEnableRollback": False,
        "SpikeEMAWindow": 50,
        "SpikeSigmaK": 10.0,
        "SpikeLogCSV": "energy/spike_trigger.csv",
        "ConfigurationLoadDir": "tpsfinal",
        "ConfigurationDumpDir": "tpsfinal",
    }

    if bc == "pbc":
        vmc_case["TRGDmin"] = d
        vmc_case["TRGDmax"] = 2 * d
        vmc_case["TRGTruncErr"] = 1e-15
    else:
        vmc_case["Dbmps_min"] = d
        vmc_case["Dbmps_max"] = 2 * d
        vmc_case["BMPSTruncErr"] = 1e-10
        vmc_case["MPSCompressScheme"] = "SVD"

    vmc = {"CaseParams": vmc_case}

    # Measure params
    meas_case: Dict[str, Any] = {
        "ThreadNum": 1,
        "InitialConfigStrategy": "Neel",
        "MC_total_samples": meas_samples,
        "WarmUp": 100,
        "MCLocalUpdateSweepsBetweenSample": 1,
        "ConfigurationLoadDir": "tpsfinal",
        "ConfigurationDumpDir": "tpsfinal",
    }

    if bc == "pbc":
        meas_case["TRGDmin"] = d
        meas_case["TRGDmax"] = 2 * d
        meas_case["TRGTruncErr"] = 1e-15
    else:
        meas_case["Dbmps_min"] = d
        meas_case["Dbmps_max"] = 2 * d
        meas_case["BMPSTruncErr"] = 1e-10
        meas_case["MPSCompressScheme"] = "SVD"

    measure = {"CaseParams": meas_case}

    # Apply overrides to vmc params
    for key, val in overrides.items():
        try:
            typed_val: Any = int(val)
        except ValueError:
            try:
                typed_val = float(val)
            except ValueError:
                typed_val = val
        vmc["CaseParams"][key] = typed_val

    return {
        "physics_unit": physics_unit,
        "physics_full": physics_full,
        "simple_update": simple_update,
        "vmc": vmc,
        "measure": measure,
    }
```

**Step 4: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: All 10 tests PASS.

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add cluster_run.py parameter computation with tests"
```

---

### Task 3: Slurm template generation

Add functions that produce slurm script content as strings, and the
`new-case` command that writes all files to disk.

**Files:**
- Modify: `scripts/cluster_run.py`
- Modify: `tests/test_cluster_run.py`

**Step 1: Write failing tests for slurm generation and new-case**

Append to `tests/test_cluster_run.py`:

```python
from cluster_run import generate_slurm_stage1, generate_slurm_stage2, generate_slurm_stage3


class TestSlurmGeneration:
    def test_stage1_contains_simple_update_and_tile(self):
        content = generate_slurm_stage1(
            job_name="hv_8x8D8_su",
            partition="256G56c",
            ntasks=1,
            walltime="24:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260217",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
            unit_lx=2, unit_ly=2, full_lx=8, full_ly=8,
        )
        assert "simple_update" in content
        assert "sitps_tile" in content
        assert "source" in content  # env loading
        assert "#SBATCH" in content

    def test_stage2_is_mpi(self):
        content = generate_slurm_stage2(
            job_name="hv_8x8D8_vmc",
            partition="256G56c",
            ntasks=56,
            walltime="72:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260217",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
        )
        assert "mpirun" in content
        assert "vmc_optimize" in content

    def test_stage3_copies_tpsfinal(self):
        content = generate_slurm_stage3(
            job_name="hv_8x8D8_meas",
            partition="256G56c",
            ntasks=56,
            walltime="72:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260217",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
        )
        assert "mc_measure" in content
        assert "cp" in content or "rsync" in content  # copies tpsfinal


class TestNewCaseFileGeneration:
    def test_new_case_creates_directory_structure(self, tmp_path):
        from cluster_run import generate_case_files

        case_dir = tmp_path / "run" / "8x8J2=0D8" / "20260217_120000"
        params = build_params(
            lx=8, ly=8, d=8, j2=0.0, bc="pbc",
            unit_lx=2, unit_ly=2,
            ntasks=56, optimizer="SR", lr=0.1,
            iterations=30,
            taus="0.5,0.2,0.1,0.05,0.02",
            step_cap=200,
            overrides={},
        )
        generate_case_files(
            case_dir=case_dir,
            params=params,
            job_name_prefix="hv_8x8D8",
            partition="256G56c",
            ntasks=56,
            walltime="72:00:00",
            remote_repo="/share/home/wanghx/HeisenbergVMCPEPS",
            unit_lx=2, unit_ly=2, full_lx=8, full_ly=8,
        )

        assert (case_dir / "params" / "physics_unit.json").exists()
        assert (case_dir / "params" / "physics_full.json").exists()
        assert (case_dir / "params" / "simple_update.json").exists()
        assert (case_dir / "params" / "vmc.json").exists()
        assert (case_dir / "params" / "measure.json").exists()
        assert (case_dir / "stage1_su_tile.slurm").exists()
        assert (case_dir / "stage2_vmc.slurm").exists()
        assert (case_dir / "stage3_measure.slurm").exists()
        assert (case_dir / "meta.json").exists()

        # Verify JSON is valid
        with open(case_dir / "params" / "vmc.json") as f:
            vmc = json.load(f)
        assert vmc["CaseParams"]["TRGDmin"] == 8
```

**Step 2: Run tests to verify they fail**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: `ImportError` for the new functions.

**Step 3: Implement slurm generators and `generate_case_files`**

Add to `scripts/cluster_run.py`:

```python
# ---------------------------------------------------------------------------
# Slurm template generation
# ---------------------------------------------------------------------------

_SLURM_HEADER = """\
#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task=1
#SBATCH --time={walltime}
#SBATCH --output={case_dir}/logs/slurm_%j_{stage}.out
#SBATCH --error={case_dir}/logs/slurm_%j_{stage}.err

set -euo pipefail
set -x

REPO={repo}
CASE_DIR={case_dir}

set +u
set +x
source /share/home/wanghx/myenv.sh
set -x
set -u

export LIBRARY_PATH="${{LD_LIBRARY_PATH:-}}"
export OMP_NUM_THREADS=1

cd "$CASE_DIR"
"""


def generate_slurm_stage1(
    job_name: str, partition: str, ntasks: int, walltime: str,
    case_dir: str, repo: str,
    unit_lx: int, unit_ly: int, full_lx: int, full_ly: int,
) -> str:
    header = _SLURM_HEADER.format(
        job_name=job_name, partition=partition, ntasks=1,
        walltime=walltime, case_dir=case_dir, repo=repo, stage="su_tile",
    )
    body = f"""
echo "=== Stage 1: Simple Update (unit cell) ==="
mkdir -p "$CASE_DIR/su_unit"
cd "$CASE_DIR/su_unit"
"$REPO/build/simple_update" "$CASE_DIR/params/physics_unit.json" "$CASE_DIR/params/simple_update.json"

echo "=== Stage 1b: Tile to full lattice ==="
"$REPO/build/sitps_tile" \\
  --input-dir "$CASE_DIR/su_unit/tpsfinal" \\
  --output-dir "$CASE_DIR/tpsfinal_tiled" \\
  --target-ly {full_ly} \\
  --target-lx {full_lx} \\
  --unit-ly {unit_ly} \\
  --unit-lx {unit_lx} \\
  --physics-json "$CASE_DIR/params/physics_full.json"

mkdir -p "$CASE_DIR/vmc"
cp -a "$CASE_DIR/tpsfinal_tiled" "$CASE_DIR/vmc/tpsfinal"

echo "[Stage 1] completed at $(date "+%F %T")"
"""
    return header + body


def generate_slurm_stage2(
    job_name: str, partition: str, ntasks: int, walltime: str,
    case_dir: str, repo: str,
) -> str:
    header = _SLURM_HEADER.format(
        job_name=job_name, partition=partition, ntasks=ntasks,
        walltime=walltime, case_dir=case_dir, repo=repo, stage="vmc",
    )
    body = """
echo "=== Stage 2: VMC Optimize ==="
cd "$CASE_DIR/vmc"
mpirun -n {ntasks} "$REPO/build/vmc_optimize" "$CASE_DIR/params/physics_full.json" "$CASE_DIR/params/vmc.json"

echo "[Stage 2] completed at $(date "+%F %T")"
""".format(ntasks=ntasks)
    return header + body


def generate_slurm_stage3(
    job_name: str, partition: str, ntasks: int, walltime: str,
    case_dir: str, repo: str,
) -> str:
    header = _SLURM_HEADER.format(
        job_name=job_name, partition=partition, ntasks=ntasks,
        walltime=walltime, case_dir=case_dir, repo=repo, stage="meas",
    )
    body = """
echo "=== Stage 3: MC Measure ==="
mkdir -p "$CASE_DIR/measure"
cp -a "$CASE_DIR/vmc/tpsfinal" "$CASE_DIR/measure/tpsfinal"
cd "$CASE_DIR/measure"
mpirun -n {ntasks} "$REPO/build/mc_measure" "$CASE_DIR/params/physics_full.json" "$CASE_DIR/params/measure.json"

echo "[Stage 3] completed at $(date "+%F %T")"
""".format(ntasks=ntasks)
    return header + body


# ---------------------------------------------------------------------------
# File generation
# ---------------------------------------------------------------------------

def generate_case_files(
    case_dir: Path,
    params: Dict[str, Any],
    job_name_prefix: str,
    partition: str,
    ntasks: int,
    walltime: str,
    remote_repo: str,
    unit_lx: int,
    unit_ly: int,
    full_lx: int,
    full_ly: int,
) -> None:
    """Write all param files, slurm scripts, and meta.json to case_dir."""

    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "params").mkdir(exist_ok=True)
    (case_dir / "logs").mkdir(exist_ok=True)
    (case_dir / "su_unit").mkdir(exist_ok=True)
    (case_dir / "vmc").mkdir(exist_ok=True)
    (case_dir / "measure").mkdir(exist_ok=True)

    remote_case = str(Path(remote_repo) / "run" / case_dir.parent.name / case_dir.name)

    # Write JSON param files
    for name, data in params.items():
        with open(case_dir / "params" / f"{name}.json", "w") as f:
            json.dump(data, f, indent=2)
            f.write("\n")

    # Write slurm scripts
    s1 = generate_slurm_stage1(
        f"{job_name_prefix}_su", partition, 1, walltime,
        remote_case, remote_repo, unit_lx, unit_ly, full_lx, full_ly,
    )
    (case_dir / "stage1_su_tile.slurm").write_text(s1)

    s2 = generate_slurm_stage2(
        f"{job_name_prefix}_vmc", partition, ntasks, walltime,
        remote_case, remote_repo,
    )
    (case_dir / "stage2_vmc.slurm").write_text(s2)

    s3 = generate_slurm_stage3(
        f"{job_name_prefix}_meas", partition, ntasks, walltime,
        remote_case, remote_repo,
    )
    (case_dir / "stage3_measure.slurm").write_text(s3)

    # Write meta.json
    meta = {
        "case": case_dir.parent.name,
        "run_id": case_dir.name,
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "partition": partition,
        "ntasks": ntasks,
        "walltime": walltime,
        "remote_repo": remote_repo,
        "params_summary": {
            "lx": full_lx, "ly": full_ly,
            "D": params["simple_update"]["CaseParams"]["Dmin"],
            "J2": params["physics_full"]["CaseParams"]["J2"],
            "bc": "pbc" if params["physics_full"]["CaseParams"]["BoundaryCondition"] == "Periodic" else "obc",
            "optimizer": params["vmc"]["CaseParams"]["OptimizerType"],
            "lr": params["vmc"]["CaseParams"]["LearningRate"],
            "iterations": params["vmc"]["CaseParams"]["MaxIterations"],
            "vmc_samples": params["vmc"]["CaseParams"]["MC_total_samples"],
            "meas_samples": params["measure"]["CaseParams"]["MC_total_samples"],
        },
        "job_ids": {},
    }
    with open(case_dir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")
```

**Step 4: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: All tests PASS.

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add slurm template generation and case file writer"
```

---

### Task 4: CLI entry points (`new-case`, `submit`, `status`, `fetch`)

Wire up the argparse CLI with all four subcommands.

**Files:**
- Modify: `scripts/cluster_run.py`

**Step 1: Add argparse main and subcommands**

Append to `scripts/cluster_run.py`:

```python
# ---------------------------------------------------------------------------
# CLI summary printer
# ---------------------------------------------------------------------------

def print_summary(case_name: str, params: Dict[str, Any], partition: str,
                  ntasks: int, walltime: str) -> None:
    su = params["simple_update"]["CaseParams"]
    vmc = params["vmc"]["CaseParams"]
    meas = params["measure"]["CaseParams"]
    pf = params["physics_full"]["CaseParams"]
    pu = params["physics_unit"]["CaseParams"]

    print(f"\n=== Case: {case_name} ===\n")
    print("Physics:")
    print(f"  Lattice:     {pf['Lx']}x{pf['Ly']} (unit cell: {pu['Lx']}x{pu['Ly']})")
    print(f"  Model:       {pf['ModelType']}")
    print(f"  J2:          {pf['J2']}")
    print(f"  BC:          {pf['BoundaryCondition']}")
    print(f"  D_peps:      {su['Dmin']}")
    print()
    print("Simple Update:")
    print(f"  Tau schedule: {su['TauScheduleTaus']}")
    print(f"  Step cap:     {su['Step']} per stage")
    print(f"  Advanced stop: enabled")
    print()
    print("VMC:")
    print(f"  Optimizer:    {vmc['OptimizerType']}")
    print(f"  LR:           {vmc['LearningRate']}")
    print(f"  Iterations:   {vmc['MaxIterations']}")
    print(f"  MC samples:   {vmc['MC_total_samples']}")
    bc_key = "TRGDmin" if "TRGDmin" in vmc else "Dbmps_min"
    bc_max_key = "TRGDmax" if "TRGDmax" in vmc else "Dbmps_max"
    print(f"  Contraction:  Dmin={vmc[bc_key]}, Dmax={vmc[bc_max_key]}")
    print()
    print("Measurement:")
    print(f"  MC samples:   {meas['MC_total_samples']}")
    print()
    print("Cluster:")
    print(f"  Partition:    {partition}")
    print(f"  Tasks:        {ntasks}")
    print(f"  Wall time:    {walltime}")


# ---------------------------------------------------------------------------
# SSH helpers
# ---------------------------------------------------------------------------

def ssh_run(cluster: str, command: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["ssh", cluster, command],
        capture_output=True, text=True, timeout=30,
    )


def scp_to_remote(local_path: str, cluster: str, remote_path: str) -> None:
    subprocess.run(
        ["scp", "-r", local_path, f"{cluster}:{remote_path}"],
        check=True, timeout=120,
    )


def scp_from_remote(cluster: str, remote_path: str, local_path: str) -> None:
    Path(local_path).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["scp", "-r", f"{cluster}:{remote_path}", local_path],
        check=True, timeout=120,
    )


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_new_case(args: argparse.Namespace) -> int:
    lx, ly = parse_lattice_spec(args.lattice)
    ulx, uly = parse_lattice_spec(args.unit_cell)

    overrides = {}
    for ov in (args.override or []):
        key, val = ov.split("=", 1)
        overrides[key] = val

    params = build_params(
        lx=lx, ly=ly, d=args.D, j2=args.J2, bc=args.bc,
        unit_lx=ulx, unit_ly=uly,
        ntasks=args.ntasks, optimizer=args.optimizer, lr=args.lr,
        iterations=args.iterations, taus=args.taus,
        step_cap=args.step_cap, overrides=overrides,
    )

    case_name = compute_case_name(lx, ly, args.J2, args.D)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    case_dir = Path("run") / case_name / timestamp

    print_summary(case_name, params, args.partition, args.ntasks, args.walltime)

    answer = input(f"\nGenerate to {case_dir}/ ? [Y/n] ").strip().lower()
    if answer and answer != "y":
        print("Aborted.")
        return 1

    job_name_prefix = f"hv_{case_name}"
    generate_case_files(
        case_dir=case_dir,
        params=params,
        job_name_prefix=job_name_prefix,
        partition=args.partition,
        ntasks=args.ntasks,
        walltime=args.walltime,
        remote_repo=args.remote_repo,
        unit_lx=ulx, unit_ly=uly,
        full_lx=lx, full_ly=ly,
    )
    print(f"\nGenerated: {case_dir}/")
    return 0


def cmd_submit(args: argparse.Namespace) -> int:
    case_dir = Path(args.case)
    if not case_dir.exists():
        print(f"ERROR: {case_dir} does not exist locally", file=sys.stderr)
        return 1

    meta_path = case_dir / "meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    remote_repo = meta["remote_repo"]
    remote_case = f"{remote_repo}/run/{meta['case']}/{meta['run_id']}"
    cluster = args.cluster

    print(f"Uploading {case_dir} to {cluster}:{remote_case} ...")
    ssh_run(cluster, f"mkdir -p {remote_case}")
    scp_to_remote(str(case_dir) + "/", cluster, remote_case)

    # Submit stage 1
    result1 = ssh_run(cluster, f"cd {remote_case} && sbatch stage1_su_tile.slurm")
    if result1.returncode != 0:
        print(f"ERROR submitting stage1: {result1.stderr}", file=sys.stderr)
        return 1
    try:
        job1 = parse_sbatch_job_id(result1.stdout)
    except ValueError as exc:
        print(f"ERROR parsing stage1 job ID: {exc}", file=sys.stderr)
        return 1
    print(f"Stage 1 (SU+tile): job {job1}")

    # Submit stage 2 with dependency
    result2 = ssh_run(cluster, f"cd {remote_case} && sbatch --dependency=afterok:{job1} stage2_vmc.slurm")
    if result2.returncode != 0:
        print(f"ERROR submitting stage2: {result2.stderr}", file=sys.stderr)
        return 1
    try:
        job2 = parse_sbatch_job_id(result2.stdout)
    except ValueError as exc:
        print(f"ERROR parsing stage2 job ID: {exc}", file=sys.stderr)
        return 1
    print(f"Stage 2 (VMC):     job {job2} (after {job1})")

    # Submit stage 3 with dependency
    result3 = ssh_run(cluster, f"cd {remote_case} && sbatch --dependency=afterok:{job2} stage3_measure.slurm")
    if result3.returncode != 0:
        print(f"ERROR submitting stage3: {result3.stderr}", file=sys.stderr)
        return 1
    try:
        job3 = parse_sbatch_job_id(result3.stdout)
    except ValueError as exc:
        print(f"ERROR parsing stage3 job ID: {exc}", file=sys.stderr)
        return 1
    print(f"Stage 3 (measure): job {job3} (after {job2})")

    # Update meta
    meta["job_ids"] = {"stage1": job1, "stage2": job2, "stage3": job3}
    meta["submitted_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")

    print(f"\nAll 3 stages submitted. Job IDs saved to {meta_path}")
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    result = ssh_run(args.cluster, f"squeue -u {args.user}")
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        return 1
    print(result.stdout)
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    case_dir = Path(args.case)
    meta_path = case_dir / "meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    remote_repo = meta["remote_repo"]
    remote_case = f"{remote_repo}/run/{meta['case']}/{meta['run_id']}"
    cluster = args.cluster
    output_dir = Path(args.output_dir)

    # Fetch energy trajectory
    energy_remote = f"{remote_case}/vmc/energy"
    energy_local = output_dir / meta["case"] / meta["run_id"] / "energy"
    print(f"Fetching energy trajectory...")
    try:
        scp_from_remote(cluster, energy_remote, str(energy_local))
        print(f"  -> {energy_local}/")
    except subprocess.CalledProcessError:
        print("  -> energy/ not found (VMC may not have run yet)")

    # Fetch measurement stats
    stats_remote = f"{remote_case}/measure/stats"
    stats_local = output_dir / meta["case"] / meta["run_id"] / "stats"
    print(f"Fetching measurement stats...")
    try:
        scp_from_remote(cluster, stats_remote, str(stats_local))
        print(f"  -> {stats_local}/")
    except subprocess.CalledProcessError:
        print("  -> stats/ not found (measurement may not have run yet)")

    # Fetch logs
    logs_remote = f"{remote_case}/logs"
    logs_local = output_dir / meta["case"] / meta["run_id"] / "logs"
    print(f"Fetching logs...")
    try:
        scp_from_remote(cluster, logs_remote, str(logs_local))
        print(f"  -> {logs_local}/")
    except subprocess.CalledProcessError:
        print("  -> logs/ not found")

    print(f"\nResults saved to {output_dir}/{meta['case']}/{meta['run_id']}/")
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="VMC-PEPS cluster run manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # new-case
    p_new = sub.add_parser("new-case", help="Generate a new case directory")
    p_new.add_argument("--lattice", required=True, help="Lattice size, e.g. 8x8")
    p_new.add_argument("--D", type=int, required=True, help="PEPS bond dimension")
    p_new.add_argument("--J2", type=float, default=0.0, help="J2 coupling (default: 0)")
    p_new.add_argument("--bc", choices=["pbc", "obc"], default="pbc", help="Boundary condition")
    p_new.add_argument("--unit-cell", default="2x2", help="Unit cell for simple update (default: 2x2)")
    p_new.add_argument("--cluster", default="susphy", help="Cluster SSH alias")
    p_new.add_argument("--partition", default="256G56c", help="Slurm partition")
    p_new.add_argument("--ntasks", type=int, default=56, help="MPI ranks")
    p_new.add_argument("--walltime", default="72:00:00", help="Wall time limit")
    p_new.add_argument("--optimizer", default="SR", choices=["SR", "Adam", "AdaGrad", "SGD"])
    p_new.add_argument("--lr", type=float, default=0.1, help="Learning rate")
    p_new.add_argument("--iterations", type=int, default=30, help="VMC iterations")
    p_new.add_argument("--taus", default="0.5,0.2,0.1,0.05,0.02", help="Tau schedule CSV")
    p_new.add_argument("--step-cap", type=int, default=200, help="Step cap per tau stage")
    p_new.add_argument("--override", action="append", help="Expert override key=value")
    p_new.add_argument("--remote-repo", default="/share/home/wanghx/HeisenbergVMCPEPS")

    # submit
    p_sub = sub.add_parser("submit", help="Upload and submit all stages")
    p_sub.add_argument("--case", required=True, help="Local case directory path")
    p_sub.add_argument("--cluster", default="susphy")

    # status
    p_st = sub.add_parser("status", help="Check running jobs")
    p_st.add_argument("--cluster", default="susphy")
    p_st.add_argument("--user", default="wanghx")

    # fetch
    p_fetch = sub.add_parser("fetch", help="Fetch results to local")
    p_fetch.add_argument("--case", required=True, help="Local case directory path")
    p_fetch.add_argument("--cluster", default="susphy")
    p_fetch.add_argument("--output-dir", default="data")

    args = parser.parse_args()

    commands = {
        "new-case": cmd_new_case,
        "submit": cmd_submit,
        "status": cmd_status,
        "fetch": cmd_fetch,
    }
    return commands[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
```

**Step 2: Smoke test the CLI**

```bash
python3 scripts/cluster_run.py --help
python3 scripts/cluster_run.py new-case --help
```

Expected: help output with all options listed.

**Step 3: Run full test suite**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: All tests PASS.

**Step 4: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add CLI entry points (new-case, submit, status, fetch)"
```

---

### Task 4b: Error path and edge case tests

Add tests for SSH/scp failures, malformed sbatch output, invalid input parsing,
and timeout handling. These use `unittest.mock` to simulate network errors
without real SSH.

**Files:**
- Modify: `tests/test_cluster_run.py`

**Step 1: Write error path tests**

Append to `tests/test_cluster_run.py`:

```python
import json
import subprocess
from unittest.mock import patch, MagicMock
from cluster_run import (
    ssh_run, cmd_submit, cmd_status, parse_sbatch_job_id,
)


class TestSSHErrorPaths:
    def test_ssh_timeout_raises(self):
        with patch("cluster_run.subprocess.run", side_effect=subprocess.TimeoutExpired("ssh", 30)):
            with pytest.raises(subprocess.TimeoutExpired):
                ssh_run("susphy", "hostname")

    def test_ssh_connection_refused(self):
        mock_result = MagicMock()
        mock_result.returncode = 255
        mock_result.stdout = ""
        mock_result.stderr = "ssh: connect to host 10.20.26.130 port 22: Connection refused"
        with patch("cluster_run.subprocess.run", return_value=mock_result):
            result = ssh_run("susphy", "hostname")
            assert result.returncode == 255


class TestSubmitErrorPaths:
    def test_submit_nonexistent_case(self, tmp_path):
        args = MagicMock()
        args.case = str(tmp_path / "nonexistent")
        args.cluster = "susphy"
        assert cmd_submit(args) == 1

    def test_submit_sbatch_parse_failure(self, tmp_path):
        """sbatch returns 0 but with unexpected output format."""
        case_dir = tmp_path / "run" / "4x4J2=0D2" / "20260217_120000"
        case_dir.mkdir(parents=True)
        meta = {
            "case": "4x4J2=0D2", "run_id": "20260217_120000",
            "remote_repo": "/share/home/wanghx/HeisenbergVMCPEPS",
            "job_ids": {},
        }
        (case_dir / "meta.json").write_text(json.dumps(meta))

        mock_mkdir = MagicMock(returncode=0, stdout="", stderr="")
        mock_bad_sbatch = MagicMock(
            returncode=0,
            stdout="sbatch: WARNING: some warning\nUnexpected output\n",
            stderr="",
        )
        with patch("cluster_run.ssh_run", side_effect=[mock_mkdir, mock_bad_sbatch]):
            with patch("cluster_run.scp_to_remote"):
                args = MagicMock()
                args.case = str(case_dir)
                args.cluster = "susphy"
                assert cmd_submit(args) == 1  # should fail gracefully


class TestStatusErrorPaths:
    def test_status_ssh_failure(self):
        mock_result = MagicMock(returncode=1, stdout="", stderr="Connection timed out")
        with patch("cluster_run.ssh_run", return_value=mock_result):
            args = MagicMock()
            args.cluster = "susphy"
            args.user = "wanghx"
            assert cmd_status(args) == 1
```

**Step 2: Run tests**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: All tests PASS (including new error path tests).

**Step 3: Commit**

```bash
git add tests/test_cluster_run.py
git commit -m "test: add error path tests for SSH, sbatch parsing, and CLI edge cases"
```

---

### Task 5: Write the skill file

Create the skill that tells Claude how to use the workflow.

**Files:**
- Create: `skills/vmcpeps-cluster-workflow.md`
- Delete: `skills/susphy_cluster_run_notes.md` (content incorporated)

**Step 1: Write the skill**

Create `skills/vmcpeps-cluster-workflow.md`:

```markdown
---
name: vmcpeps-cluster-workflow
description: "Use when the user asks to run calculations, submit jobs, check
job status, fetch results, or manage cases on the susphy HPC cluster. Covers
the full pipeline: simple update, tile, VMC optimize, MC measure."
---

# VMC-PEPS Cluster Workflow

Use this skill when setting up, submitting, monitoring, or fetching results
for VMC-PEPS calculations on the susphy cluster.

## Pipeline Overview

Three slurm jobs with dependency chaining:

1. **Stage 1 (SU + tile):** Single-node. Runs `simple_update` on unit cell
   with tau-schedule + advanced-stop, then `sitps_tile` to expand to full
   lattice.
2. **Stage 2 (VMC):** MPI job. Runs `vmc_optimize` on full lattice.
   `--dependency=afterok:stage1`.
3. **Stage 3 (Measure):** MPI job. Runs `mc_measure` on optimized state.
   `--dependency=afterok:stage2`.

## CLI Reference

All commands via `python3 scripts/cluster_run.py`:

### Generate a new case

```bash
python3 scripts/cluster_run.py new-case \
  --lattice 8x8 --D 8 --J2 0.0 --bc pbc \
  --unit-cell 2x2 \
  --partition 256G56c --ntasks 56
```

Key optional flags: `--optimizer SR`, `--lr 0.1`, `--iterations 30`,
`--taus "0.5,0.2,0.1,0.05,0.02"`, `--override key=value`.

### Submit to cluster

```bash
python3 scripts/cluster_run.py submit \
  --case run/8x8J2=0D8/<timestamp>
```

### Check status

```bash
python3 scripts/cluster_run.py status
```

Or directly: `ssh susphy "squeue -u wanghx"`

### Fetch results

```bash
python3 scripts/cluster_run.py fetch \
  --case run/8x8J2=0D8/<timestamp> --output-dir data/
```

## Pre-Submit Checklist

Before submitting any job:

1. Show the full parameter summary to the user
2. Confirm the user approves the parameters (especially D, lattice size,
   optimizer, iteration count)
3. Verify SSH connectivity: `ssh susphy "hostname"`
4. Check no conflicting jobs: `ssh susphy "squeue -u wanghx"`

## Parameter Defaults

- MC samples (VMC): `max(10000, min(20000, Lx*Ly*200))`, rounded to ntasks
- MC samples (measure): 2x VMC samples
- TRGDmin = D, TRGDmax = 2*D (PBC only)
- Tau schedule: `0.5, 0.2, 0.1, 0.05, 0.02`
- Advanced stop: energy_tol=1e-6, lambda_tol=1e-5, patience=3

## Cluster Operational Notes (susphy)

### Partition selection
- Use `256G56c` for production runs. The 24c partitions have ABI
  compatibility issues with login-node builds.

### Environment loading
Slurm scripts must source the environment with nounset temporarily disabled:

```bash
set +u; set +x
source /share/home/wanghx/myenv.sh
set -x; set -u
```

### Failure signatures
- `FAILED 0:00 with empty output on 64G24c`: partition unusable, switch to
  256G56c.
- `libmkl_* not found` / `GLIBCXX_* not found`: env not loaded in job script.
- `BASHRCSOURCED unbound`: don't source .bashrc under set -u.

### Expected success checkpoints
- `Simple Update completed.`
- `SplitIndexTPS saved to: tpsfinal`
- `Optimization completed!`
- `Rank 0: statistic data finished.`

## Result Files

After a successful run, key output files are:

- `vmc/energy/energy_trajectory.csv` -- columns: iteration, energy,
  energy_error, gradient_norm
- `measure/stats/energy.csv` -- final energy: index, mean, stderr
- `measure/stats/SzSz_all2all.csv` -- spin-spin correlations
- `measure/stats/bond_energy_{h,v}_{mean,stderr}.csv` -- bond energies

## Directory Structure

```
run/<case>/<timestamp>/
  params/         -- frozen input JSONs
  stage1_su_tile.slurm
  stage2_vmc.slurm
  stage3_measure.slurm
  su_unit/        -- simple update working dir
  tpsfinal_tiled/ -- tiled TPS
  vmc/            -- VMC working dir (energy/, tpsfinal/, checkpoints/)
  measure/        -- measurement working dir (stats/, tpsfinal/)
  logs/           -- slurm stdout/stderr
  meta.json       -- run metadata + job IDs
```
```

**Step 2: Delete the old notes file**

```bash
rm skills/susphy_cluster_run_notes.md
```

**Step 3: Commit**

```bash
git add skills/vmcpeps-cluster-workflow.md
git rm skills/susphy_cluster_run_notes.md 2>/dev/null || true
git commit -m "feat: add vmcpeps-cluster-workflow skill, remove old notes"
```

---

### Task 6: End-to-end local test

Verify the CLI generates correct files by running `new-case` locally and
inspecting the output.

**Files:** None modified (integration test)

**Step 1: Generate a test case**

```bash
cd /Users/wanghaoxin/GitHub/HeisenbergVMCPEPS
python3 scripts/cluster_run.py new-case \
  --lattice 8x8 --D 8 --J2 0.0 --bc pbc \
  --unit-cell 2x2 \
  --partition 256G56c --ntasks 56
```

Type `Y` when prompted.

**Step 2: Verify generated files**

```bash
ls run/8x8J2=0D8/*/
cat run/8x8J2=0D8/*/params/simple_update.json
cat run/8x8J2=0D8/*/params/vmc.json
cat run/8x8J2=0D8/*/meta.json
cat run/8x8J2=0D8/*/stage1_su_tile.slurm
```

Verify:
- `simple_update.json` has `Dmin=8, Dmax=8, TauScheduleEnabled=true`
- `vmc.json` has `TRGDmin=8, TRGDmax=16, OptimizerType=SR`
- `meta.json` has correct case name and params summary
- `stage1_su_tile.slurm` has `simple_update` and `sitps_tile` commands

**Step 3: Clean up test output**

```bash
rm -rf run/8x8J2=0D8/
```

**Step 4: Run all tests one final time**

```bash
python3 -m pytest tests/test_cluster_run.py -v
```

Expected: All tests PASS.

**Step 5: Commit everything**

If any remaining uncommitted changes exist (e.g., test adjustments):

```bash
git add -A scripts/ tests/ skills/ docs/plans/
git commit -m "feat: cluster workflow CLI and skill - complete implementation"
```

---

### Task 7: Live cluster smoke test (optional)

Test the full submit flow on susphy with a small OBC case (OBC is faster than
PBC since it uses BMPS contraction instead of TRG).

**Step 1: Generate a small test case**

```bash
python3 scripts/cluster_run.py new-case \
  --lattice 4x4 --D 2 --J2 0.0 --bc obc \
  --unit-cell 4x4 \
  --partition 256G56c --ntasks 4 \
  --iterations 3 --step-cap 8 \
  --walltime 00:15:00
```

Note: OBC with `--unit-cell 4x4` means simple update runs on the full lattice
directly (no tiling step needed, but `sitps_tile` will be a no-op identity
copy). This is the fastest possible smoke test.

**Step 2: Submit**

```bash
python3 scripts/cluster_run.py submit \
  --case run/4x4J2=0D2/<timestamp>
```

**Step 3: Monitor**

```bash
python3 scripts/cluster_run.py status
```

Wait for all 3 stages to complete (should take under 5 minutes for D=2 4x4 OBC).

**Step 4: Fetch results**

```bash
python3 scripts/cluster_run.py fetch \
  --case run/4x4J2=0D2/<timestamp> --output-dir data/
ls data/4x4J2=0D2/*/energy/
ls data/4x4J2=0D2/*/stats/
```

**Step 5: Clean up test run**

```bash
rm -rf run/4x4J2=0D2/
ssh susphy "rm -rf ~/HeisenbergVMCPEPS/run/4x4J2=0D2/"
```
