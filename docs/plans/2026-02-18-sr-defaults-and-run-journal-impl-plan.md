# SR Defaults & Run Journal Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add battle-tested SR/CG defaults with auto-scaled WarmUp, QoS support, 2-month walltime default, and a persistent run journal that the AI consults for parameter suggestions.

**Architecture:** Pure-function additions to `scripts/cluster_run.py` (parameter computation, journal parsing), slurm template updates (QoS line), skill file enhancements (journal protocol). No new files except `run_journal.jsonl` (generated at runtime).

**Tech Stack:** Python 3 stdlib (json, csv, re, pathlib), pytest

---

### Task 1: Auto-scaled WarmUp and CGResidueRestart

**Files:**
- Modify: `scripts/cluster_run.py:44-48` (add new functions after `compute_mc_samples`)
- Modify: `scripts/cluster_run.py:129` (use `compute_warmup` in `build_params`)
- Modify: `scripts/cluster_run.py:136` (use `compute_cg_residue_restart` in `build_params`)
- Test: `tests/test_cluster_run.py`

**Step 1: Write the failing tests**

Add to `tests/test_cluster_run.py`, after the import block (line 17), add `compute_warmup, compute_cg_residue_restart` to the import. Then add test classes after `TestComputeCaseName` (line 83):

```python
from cluster_run import (
    compute_mc_samples, compute_case_name, build_params,
    parse_sbatch_job_id, parse_lattice_spec,
    compute_warmup, compute_cg_residue_restart,
    generate_slurm_stage1, generate_slurm_stage2, generate_slurm_stage3,
    generate_case_files, ssh_run, cmd_submit, cmd_status,
)

# ... existing tests ...

class TestComputeWarmup:
    def test_small_lattice_clamps_to_minimum(self):
        assert compute_warmup(4, 4) == 100  # 4*4*2=32, clamped to 100

    def test_medium_lattice_scales(self):
        assert compute_warmup(8, 8) == 128  # 8*8*2=128

    def test_large_lattice_clamps_to_maximum(self):
        assert compute_warmup(16, 16) == 500  # 16*16*2=512, clamped to 500

    def test_rectangular(self):
        assert compute_warmup(12, 8) == 192  # 12*8*2=192


class TestComputeCgResidueRestart:
    def test_small_lattice(self):
        assert compute_cg_residue_restart(8, 8) == 10

    def test_large_lattice(self):
        assert compute_cg_residue_restart(16, 16) == 20

    def test_boundary(self):
        assert compute_cg_residue_restart(16, 8) == 20  # max(16,8)>=16
        assert compute_cg_residue_restart(15, 15) == 10  # max(15,15)<16
```

**Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_cluster_run.py::TestComputeWarmup -v`
Expected: ImportError — `compute_warmup` not found

**Step 3: Write minimal implementation**

In `scripts/cluster_run.py`, after `compute_mc_samples` (line 48), add:

```python
def compute_warmup(lx: int, ly: int) -> int:
    """Auto-scale WarmUp: 2 * Lx * Ly, clamped to [100, 500]."""
    return max(100, min(500, lx * ly * 2))


def compute_cg_residue_restart(lx: int, ly: int) -> int:
    """CGResidueRestart: 10 for most systems, 20 for large (max side >= 16)."""
    return 20 if max(lx, ly) >= 16 else 10
```

Then update `build_params` — change line 129 from `"WarmUp": 100,` to:

```python
        "WarmUp": compute_warmup(lx, ly),
```

And change line 136 from `"CGResidueRestart": 10,` to:

```python
        "CGResidueRestart": compute_cg_residue_restart(lx, ly),
```

Also update the measure params — change line 178 from `"WarmUp": 100,` to:

```python
        "WarmUp": compute_warmup(lx, ly),
```

**Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_cluster_run.py -v`
Expected: All pass (existing + new)

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add auto-scaled WarmUp and CGResidueRestart for SR defaults"
```

---

### Task 2: QoS support and walltime default

**Files:**
- Modify: `scripts/cluster_run.py:219-246` (`_SLURM_HEADER` template)
- Modify: `scripts/cluster_run.py:249-254` (`generate_slurm_stage1` signature)
- Modify: `scripts/cluster_run.py:289-292` (`generate_slurm_stage2` signature)
- Modify: `scripts/cluster_run.py:307-310` (`generate_slurm_stage3` signature)
- Modify: `scripts/cluster_run.py:331-343` (`generate_case_files` signature)
- Modify: `scripts/cluster_run.py:660` (walltime default)
- Modify: `scripts/cluster_run.py:650-667` (add `--qos` arg)
- Test: `tests/test_cluster_run.py`

**Step 1: Write the failing tests**

Add to test file:

```python
class TestQosInSlurm:
    def test_qos_appears_in_stage1(self):
        content = generate_slurm_stage1(
            job_name="hv_8x8D8_su", partition="256G56c", ntasks=1,
            walltime="60-00:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260218",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
            unit_lx=2, unit_ly=2, full_lx=8, full_ly=8,
            bc="pbc", qos="wanghxqos",
        )
        assert "#SBATCH --qos=wanghxqos" in content

    def test_qos_appears_in_stage2(self):
        content = generate_slurm_stage2(
            job_name="hv_8x8D8_vmc", partition="256G56c", ntasks=56,
            walltime="60-00:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260218",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
            qos="wanghxqos",
        )
        assert "#SBATCH --qos=wanghxqos" in content

    def test_no_qos_when_empty(self):
        content = generate_slurm_stage2(
            job_name="hv_8x8D8_vmc", partition="256G56c", ntasks=56,
            walltime="60-00:00:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/8x8J2=0D8/20260218",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
        )
        assert "--qos" not in content
```

**Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_cluster_run.py::TestQosInSlurm -v`
Expected: TypeError — unexpected keyword argument `qos`

**Step 3: Write minimal implementation**

Add `qos` parameter to `_SLURM_HEADER`. Change the template to accept an optional QoS line. The simplest approach: add `qos` param to the header format function.

Replace `_SLURM_HEADER` (lines 219-246) with:

```python
def _slurm_header(
    job_name: str, partition: str, ntasks: int, walltime: str,
    case_dir: str, repo: str, stage: str, qos: str = "",
) -> str:
    qos_line = f"\n#SBATCH --qos={qos}" if qos else ""
    return f"""\
#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task=1
#SBATCH --time={walltime}{qos_line}
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
```

Then update all three `generate_slurm_stage*` functions to accept `qos: str = ""` and pass it through to `_slurm_header()`. Update `generate_case_files` to accept and pass `qos`. Update all call sites.

Change `--walltime` default (line 660) from `"72:00:00"` to `"60-00:00:00"`.

Add `--qos` argument after `--walltime`:

```python
p_new.add_argument("--qos", default="wanghxqos", help="Slurm QoS (default: wanghxqos)")
```

Update existing slurm generation tests that call these functions to still work (they don't pass `qos`, so the default `""` keeps them passing).

**Step 4: Run tests**

Run: `python3 -m pytest tests/test_cluster_run.py -v`
Expected: All pass

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add QoS support and change walltime default to 60 days"
```

---

### Task 3: WarmUp hint in print_summary

**Files:**
- Modify: `scripts/cluster_run.py:410-446` (`print_summary`)

**Step 1: Update `print_summary`**

This is a display-only change. In `print_summary`, add WarmUp and CG info lines. After the `MC samples` line (line 435), add:

```python
    print(f"  WarmUp:       {vmc['WarmUp']} (auto: 2 x Lx x Ly, range [100, 500])")
    print(f"  CG:           MaxIter={vmc['CGMaxIter']}, Tol={vmc['CGTol']}, "
          f"DiagShift={vmc['CGDiagShift']}, Restart={vmc['CGResidueRestart']}")
```

Also update the Cluster section to show QoS. After `Wall time` line (446), add:

```python
    if hasattr(args_ref, 'qos') if args_ref else False:
        print(f"  QoS:          {qos}")
```

Actually, simpler: add `qos` as a parameter to `print_summary`:

```python
def print_summary(case_name: str, params: Dict[str, Any], partition: str,
                  ntasks: int, walltime: str, qos: str = "") -> None:
```

And at the end of the Cluster section:

```python
    if qos:
        print(f"  QoS:          {qos}")
```

**Step 2: No test needed** — display-only, verified by existing E2E test flow.

**Step 3: Commit**

```bash
git add scripts/cluster_run.py
git commit -m "feat: show WarmUp hint and CG defaults in case summary"
```

---

### Task 4: Journal entry builder and parser

**Files:**
- Modify: `scripts/cluster_run.py` (add new functions after SSH helpers section, ~line 472)
- Test: `tests/test_cluster_run.py`

**Step 1: Write the failing tests**

```python
from cluster_run import (
    # ... existing imports ...,
    compute_warmup, compute_cg_residue_restart,
    parse_energy_trajectory, build_journal_entry, read_journal,
    suggest_walltime,
    # ...
)


class TestParseEnergyTrajectory:
    def test_parses_csv(self, tmp_path):
        csv = tmp_path / "energy_trajectory.csv"
        csv.write_text(
            "iteration,energy,energy_err,gradient_norm\n"
            "0,-8.5,0.02,0.44\n"
            "1,-8.56,0.018,0.23\n"
            "2,-8.62,0.017,0.17\n"
        )
        result = parse_energy_trajectory(str(csv))
        assert result["final_energy"] == pytest.approx(-8.62)
        assert result["best_energy"] == pytest.approx(-8.62)
        assert result["best_iter"] == 2
        assert result["total_iterations"] == 3
        assert result["gradient_norm_final"] == pytest.approx(0.17)

    def test_missing_file_returns_none(self):
        assert parse_energy_trajectory("/nonexistent/path.csv") is None


class TestBuildJournalEntry:
    def test_builds_entry_from_meta_and_energy(self, tmp_path):
        meta = {
            "case": "4x4J2=0D2", "run_id": "20260218_120000",
            "partition": "256G56c", "ntasks": 4,
            "walltime": "60-00:00:00",
            "remote_repo": "/share/home/wanghx/HeisenbergVMCPEPS",
            "params_summary": {
                "lx": 4, "ly": 4, "D": 2, "J2": 0.0, "bc": "obc",
                "optimizer": "SR", "lr": 0.1, "iterations": 3,
                "vmc_samples": 10000, "meas_samples": 20000,
            },
        }
        energy_data = {
            "final_energy": -8.62, "final_energy_err": 0.017,
            "best_energy": -8.62, "best_iter": 2,
            "gradient_norm_final": 0.17, "spike_count": 0,
            "total_iterations": 3,
        }
        entry = build_journal_entry(meta, energy_data, walltime_used_seconds=52)
        assert entry["case"] == "4x4J2=0D2"
        assert entry["outcome"]["final_energy"] == -8.62
        assert entry["hpc"]["walltime_used_seconds"] == 52
        assert entry["insights"] == []


class TestReadJournal:
    def test_reads_jsonl(self, tmp_path):
        jf = tmp_path / "run_journal.jsonl"
        jf.write_text(
            '{"case":"4x4J2=0D2","outcome":{"final_energy":-8.62}}\n'
            '{"case":"8x8J2=0D8","outcome":{"final_energy":-25.1}}\n'
        )
        entries = read_journal(str(jf))
        assert len(entries) == 2
        assert entries[1]["case"] == "8x8J2=0D8"

    def test_missing_file_returns_empty(self, tmp_path):
        entries = read_journal(str(tmp_path / "nonexistent.jsonl"))
        assert entries == []


class TestSuggestWalltime:
    def test_no_history_returns_none(self):
        assert suggest_walltime([], 8, 8, 8, "pbc") is None

    def test_with_similar_runs(self):
        journal = [
            {"params": {"lx": 8, "ly": 8, "D": 8, "bc": "pbc"},
             "hpc": {"walltime_used_seconds": 1000}},
            {"params": {"lx": 8, "ly": 8, "D": 8, "bc": "pbc"},
             "hpc": {"walltime_used_seconds": 2000}},
            {"params": {"lx": 8, "ly": 8, "D": 8, "bc": "pbc"},
             "hpc": {"walltime_used_seconds": 3000}},
        ]
        result = suggest_walltime(journal, 8, 8, 8, "pbc")
        # 5x median(1000,2000,3000) = 5*2000 = 10000 seconds
        assert result == 10000

    def test_ignores_different_config(self):
        journal = [
            {"params": {"lx": 4, "ly": 4, "D": 2, "bc": "obc"},
             "hpc": {"walltime_used_seconds": 100}},
        ]
        assert suggest_walltime(journal, 8, 8, 8, "pbc") is None
```

**Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_cluster_run.py::TestParseEnergyTrajectory -v`
Expected: ImportError

**Step 3: Write minimal implementation**

Add after the SSH helpers section in `cluster_run.py` (~line 472):

```python
# ---------------------------------------------------------------------------
# Run journal
# ---------------------------------------------------------------------------

import csv
import statistics

def parse_energy_trajectory(csv_path: str) -> dict | None:
    """Parse energy_trajectory.csv into summary dict. Returns None if file missing."""
    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
    except (FileNotFoundError, OSError):
        return None
    if not rows:
        return None

    energies = [(int(r["iteration"]), float(r["energy"])) for r in rows]
    energy_errs = [float(r["energy_err"]) for r in rows]
    grad_norms = [float(r["gradient_norm"]) for r in rows]

    best_iter, best_energy = min(energies, key=lambda x: x[1])
    final_energy = energies[-1][1]
    final_energy_err = energy_errs[-1]

    # Count spikes: energy increase > 10x median |delta|
    deltas = [abs(energies[i+1][1] - energies[i][1]) for i in range(len(energies)-1)]
    spike_count = 0
    if deltas:
        med_delta = statistics.median(deltas) if deltas else 0
        threshold = max(med_delta * 10, 1e-6)
        spike_count = sum(1 for d in deltas if d > threshold)

    return {
        "final_energy": final_energy,
        "final_energy_err": final_energy_err,
        "best_energy": best_energy,
        "best_iter": best_iter,
        "gradient_norm_final": grad_norms[-1],
        "spike_count": spike_count,
        "total_iterations": len(rows),
    }


def build_journal_entry(
    meta: dict,
    energy_data: dict | None,
    walltime_used_seconds: int | None = None,
) -> dict:
    """Build a journal entry from meta.json and parsed energy data."""
    ps = meta.get("params_summary", {})
    entry = {
        "case": meta["case"],
        "run_id": meta["run_id"],
        "timestamp": meta.get("created_at", ""),
        "params": {
            "lx": ps.get("lx"), "ly": ps.get("ly"),
            "D": ps.get("D"), "J2": ps.get("J2"),
            "bc": ps.get("bc"),
            "optimizer": ps.get("optimizer"),
            "lr": ps.get("lr"),
            "iterations": ps.get("iterations"),
            "vmc_samples": ps.get("vmc_samples"),
        },
        "outcome": energy_data or {},
        "hpc": {
            "partition": meta.get("partition", ""),
            "ntasks": meta.get("ntasks", 0),
            "walltime_requested": meta.get("walltime", ""),
            "walltime_used_seconds": walltime_used_seconds,
        },
        "insights": [],
    }
    return entry


def read_journal(path: str = "run_journal.jsonl") -> list:
    """Read all entries from the run journal. Returns [] if file missing."""
    try:
        with open(path) as f:
            return [json.loads(line) for line in f if line.strip()]
    except (FileNotFoundError, OSError):
        return []


def append_journal(entry: dict, path: str = "run_journal.jsonl") -> None:
    """Append a single entry to the run journal."""
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def suggest_walltime(
    journal: list, lx: int, ly: int, d: int, bc: str,
) -> int | None:
    """Suggest walltime (seconds) as 5x median of similar past runs.
    Returns None if no similar runs found."""
    similar = [
        e["hpc"]["walltime_used_seconds"]
        for e in journal
        if (e.get("params", {}).get("lx") == lx
            and e.get("params", {}).get("ly") == ly
            and e.get("params", {}).get("D") == d
            and e.get("params", {}).get("bc") == bc
            and e.get("hpc", {}).get("walltime_used_seconds") is not None)
    ]
    if not similar:
        return None
    return int(5 * statistics.median(similar))
```

**Step 4: Run tests**

Run: `python3 -m pytest tests/test_cluster_run.py -v`
Expected: All pass

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: add journal entry builder, energy parser, and walltime suggestion"
```

---

### Task 5: Integrate journal into fetch command

**Files:**
- Modify: `scripts/cluster_run.py:597-636` (`cmd_fetch`)

**Step 1: Write the failing test**

```python
class TestFetchAppendsJournal:
    def test_fetch_creates_journal_entry(self, tmp_path):
        # Set up case directory with meta.json
        case_dir = tmp_path / "run" / "4x4J2=0D2" / "20260218_120000"
        case_dir.mkdir(parents=True)
        meta = {
            "case": "4x4J2=0D2", "run_id": "20260218_120000",
            "partition": "256G56c", "ntasks": 4,
            "walltime": "60-00:00:00",
            "remote_repo": "/share/home/wanghx/HeisenbergVMCPEPS",
            "params_summary": {
                "lx": 4, "ly": 4, "D": 2, "J2": 0.0, "bc": "obc",
                "optimizer": "SR", "lr": 0.1, "iterations": 3,
                "vmc_samples": 10000, "meas_samples": 20000,
            },
        }
        (case_dir / "meta.json").write_text(json.dumps(meta))

        # Set up fetched energy data
        energy_dir = tmp_path / "data" / "4x4J2=0D2" / "20260218_120000" / "energy"
        energy_dir.mkdir(parents=True)
        (energy_dir / "energy_trajectory.csv").write_text(
            "iteration,energy,energy_err,gradient_norm\n"
            "0,-8.5,0.02,0.44\n"
            "1,-8.56,0.018,0.23\n"
            "2,-8.62,0.017,0.17\n"
        )

        journal_path = str(tmp_path / "run_journal.jsonl")

        # Mock scp to copy from "remote" (just use local paths)
        def mock_scp(cluster, remote, local):
            # Simulate fetching energy dir
            if "energy" in remote:
                import shutil
                Path(local).parent.mkdir(parents=True, exist_ok=True)
                shutil.copytree(str(energy_dir), local, dirs_exist_ok=True)
            else:
                raise subprocess.CalledProcessError(1, "scp")

        with patch("cluster_run.scp_from_remote", side_effect=mock_scp):
            args = MagicMock()
            args.case = str(case_dir)
            args.cluster = "susphy"
            args.output_dir = str(tmp_path / "data")
            args.journal = journal_path
            cmd_fetch(args)

        entries = read_journal(journal_path)
        assert len(entries) == 1
        assert entries[0]["case"] == "4x4J2=0D2"
        assert entries[0]["outcome"]["final_energy"] == pytest.approx(-8.62)
```

**Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_cluster_run.py::TestFetchAppendsJournal -v`
Expected: FAIL (cmd_fetch doesn't accept `args.journal` yet)

**Step 3: Update `cmd_fetch`**

At the end of `cmd_fetch`, after fetching logs, add journal append logic:

```python
    # Append to run journal
    journal_path = getattr(args, 'journal', 'run_journal.jsonl')
    energy_csv = output_dir / meta["case"] / meta["run_id"] / "energy" / "energy_trajectory.csv"
    energy_data = parse_energy_trajectory(str(energy_csv))
    entry = build_journal_entry(meta, energy_data)
    append_journal(entry, journal_path)
    if energy_data:
        print(f"\nJournal: appended entry for {meta['case']}/{meta['run_id']}")
        print(f"  Energy: {energy_data['final_energy']:.6f} ± {energy_data['final_energy_err']:.4f}")
        print(f"  Best:   {energy_data['best_energy']:.6f} at iter {energy_data['best_iter']}")
        print(f"  Spikes: {energy_data['spike_count']}")
    else:
        print(f"\nJournal: appended entry (no energy data yet)")
```

Add `--journal` argument to the fetch subparser:

```python
p_fetch.add_argument("--journal", default="run_journal.jsonl", help="Journal file path")
```

**Step 4: Run tests**

Run: `python3 -m pytest tests/test_cluster_run.py -v`
Expected: All pass

**Step 5: Commit**

```bash
git add scripts/cluster_run.py tests/test_cluster_run.py
git commit -m "feat: auto-append journal entry on fetch with energy trajectory parsing"
```

---

### Task 6: Add run_journal.jsonl to .gitignore

**Files:**
- Modify: `.gitignore:73` (add after `.worktrees/`)

**Step 1: Add the line**

After line 73 (`.worktrees/`), add:

```
run_journal.jsonl
```

**Step 2: Commit**

```bash
git add .gitignore
git commit -m "chore: gitignore run_journal.jsonl"
```

---

### Task 7: Update skill file with journal protocol

**Files:**
- Modify: `skills/vmcpeps-cluster-workflow.md`

**Step 1: Add three new sections**

After the existing "Directory Structure" section (end of file), add:

```markdown
## Run Journal Protocol

At the start of any run-related conversation:

1. Read `run_journal.jsonl` if it exists
2. Summarize relevant past runs for the current discussion
3. When proposing parameters, cite relevant history:
   "Based on run 20260217 (8x8 D=8 PBC), SR lr=0.1 converged in 22
   iters with 1 spike. Suggesting lr=0.05."

After fetching results (`fetch` command auto-appends to journal):

4. Review the new journal entry
5. Analyze convergence quality (see checklist below)
6. Suggest improvements if warranted

## Post-Fetch Analysis Checklist

After fetching results, evaluate:

- **Energy trend**: converging / plateaued / diverging?
- **Spike count**: 0 ideal, >2 suggests lr or CGDiagShift issue
- **Accept rate**: 0.2-0.3 is healthy for spin-1/2
- **Walltime utilization**: <10% means walltime too generous, suggest 5x median
- **CG iterations**: hitting CGMaxIter often? suggest increasing it
- **Gradient norm**: should decrease; flat means optimization stuck
- **EvalT / UpdateT**: compare to similar past runs, flag regressions

## Improvement Suggestion Categories

When patterns emerge across multiple runs, suggest:

**Parameter tuning:**
- "D=8 runs consistently need lr<=0.05 to avoid spikes"
- "WarmUp=128 insufficient for 12x12; energy unstable in first iterations"

**HPC resource optimization:**
- "256G56c uses only 30GB for D=5; smaller partition available?"
- "Last 3 runs used <7% walltime; suggest 5x median = 2d12h"
- "Consider --ntasks 28 with OMP_NUM_THREADS=2 for memory-bound D=10"

**Code improvements:**
- "SpikeFactorErr=100 triggers in 4/5 D>=8 runs; consider tightening default
  in enhanced_params_parser.h"
- "CGDiagShift=1e-4 causes CG stalling at D=10; scale with D?"

**Performance:**
- "EvalT increased 3x from D=5 to D=8 but UpdateT stayed flat;
  contraction is the bottleneck, not the CG solver"

## Walltime Strategy

- **No history**: Default 60 days (safe ceiling for unknown jobs)
- **With history**: AI reads journal, finds similar runs (same lattice/D/BC),
  suggests 5x median walltime. Displayed in summary.
- **User override**: `--walltime` always wins
```

**Step 2: Commit**

```bash
git add skills/vmcpeps-cluster-workflow.md
git commit -m "feat: add journal protocol, analysis checklist, and improvement categories to skill"
```

---

## Task Dependency Summary

| Task | Description | Depends On |
|------|-------------|------------|
| 1 | Auto-scaled WarmUp + CGResidueRestart | — |
| 2 | QoS support + walltime default | — |
| 3 | WarmUp hint in print_summary | 1, 2 |
| 4 | Journal entry builder + parser | — |
| 5 | Integrate journal into fetch | 4 |
| 6 | Gitignore run_journal.jsonl | — |
| 7 | Update skill file | — |

**Parallel batches:**
- Batch 1: Tasks 1, 2, 4, 6, 7 (all independent)
- Batch 2: Tasks 3, 5 (depend on Batch 1)
