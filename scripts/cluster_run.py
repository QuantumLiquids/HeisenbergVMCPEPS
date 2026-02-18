#!/usr/bin/env python3
"""CLI for generating, submitting, and managing VMC-PEPS cluster runs."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_sbatch_job_id(stdout: str) -> str:
    """Extract job ID from sbatch output using regex. Raises ValueError on failure."""
    match = re.search(r"Submitted batch job (\d+)", stdout)
    if not match:
        raise ValueError(f"Could not parse sbatch job ID from output: {stdout!r}")
    return match.group(1)


def parse_lattice_spec(spec: str) -> tuple:
    """Parse 'LxXLy' string (e.g. '8x8') into (lx, ly). Raises ValueError on bad input."""
    match = re.fullmatch(r"(\d+)[xX](\d+)", spec.strip())
    if not match:
        raise ValueError(f"Invalid lattice spec: {spec!r}. Expected format: 8x8")
    lx, ly = int(match.group(1)), int(match.group(2))
    if lx <= 0 or ly <= 0:
        raise ValueError(f"Lattice dimensions must be positive, got {lx}x{ly}")
    return lx, ly


# ---------------------------------------------------------------------------
# Parameter computation (pure functions, no I/O)
# ---------------------------------------------------------------------------

def compute_mc_samples(lx: int, ly: int, ntasks: int) -> int:
    """Compute MC_total_samples scaled by lattice size, clamped to [10000, 20000]."""
    raw = lx * ly * 200
    clamped = max(10000, min(20000, raw))
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
    bc: str = "pbc",
) -> str:
    header = _SLURM_HEADER.format(
        job_name=job_name, partition=partition, ntasks=1,
        walltime=walltime, case_dir=case_dir, repo=repo, stage="su_tile",
    )
    # sitps_tile rejects --unit-ly/--unit-lx for OBC mode
    if bc == "pbc":
        unit_args = f"""\
  --unit-ly {unit_ly} \\
  --unit-lx {unit_lx} \\
"""
    else:
        unit_args = ""
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
{unit_args}  --physics-json "$CASE_DIR/params/physics_full.json"

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

    for name, data in params.items():
        with open(case_dir / "params" / f"{name}.json", "w") as f:
            json.dump(data, f, indent=2)
            f.write("\n")

    bc = "pbc" if params["physics_full"]["CaseParams"]["BoundaryCondition"] == "Periodic" else "obc"
    s1 = generate_slurm_stage1(
        f"{job_name_prefix}_su", partition, 1, walltime,
        remote_case, remote_repo, unit_lx, unit_ly, full_lx, full_ly,
        bc=bc,
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

    energy_remote = f"{remote_case}/vmc/energy"
    energy_local = output_dir / meta["case"] / meta["run_id"] / "energy"
    print("Fetching energy trajectory...")
    try:
        scp_from_remote(cluster, energy_remote, str(energy_local))
        print(f"  -> {energy_local}/")
    except subprocess.CalledProcessError:
        print("  -> energy/ not found (VMC may not have run yet)")

    stats_remote = f"{remote_case}/measure/stats"
    stats_local = output_dir / meta["case"] / meta["run_id"] / "stats"
    print("Fetching measurement stats...")
    try:
        scp_from_remote(cluster, stats_remote, str(stats_local))
        print(f"  -> {stats_local}/")
    except subprocess.CalledProcessError:
        print("  -> stats/ not found (measurement may not have run yet)")

    logs_remote = f"{remote_case}/logs"
    logs_local = output_dir / meta["case"] / meta["run_id"] / "logs"
    print("Fetching logs...")
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
