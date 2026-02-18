"""Tests for cluster_run.py parameter computation and CLI helpers."""

import json
import subprocess
import sys
import os
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from cluster_run import (
    compute_mc_samples, compute_case_name, build_params,
    compute_warmup, compute_cg_residue_restart,
    parse_sbatch_job_id, parse_lattice_spec,
    generate_slurm_stage1, generate_slurm_stage2, generate_slurm_stage3,
    generate_case_files, ssh_run, cmd_submit, cmd_status, cmd_fetch,
    parse_energy_trajectory, build_journal_entry, read_journal, suggest_walltime,
)


# ---- Parsing helpers ----

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


# ---- Parameter computation ----

class TestComputeMcSamples:
    def test_small_lattice_clamps_to_minimum(self):
        result = compute_mc_samples(4, 4, 56)
        assert result >= 10000
        assert result % 56 == 0

    def test_medium_lattice_scales(self):
        result = compute_mc_samples(8, 8, 56)
        assert 10000 <= result <= 20000
        assert result % 56 == 0

    def test_large_lattice_clamps_to_maximum(self):
        result = compute_mc_samples(16, 16, 56)
        # After clamping to 20000, rounds up to next multiple of ntasks
        assert result <= 20000 + 56
        assert result % 56 == 0


class TestComputeCaseName:
    def test_integer_j2(self):
        assert compute_case_name(8, 8, 0.0, 8) == "8x8J2=0D8"

    def test_fractional_j2(self):
        assert compute_case_name(8, 8, 0.5, 8) == "8x8J2=0.5D8"

    def test_different_lattice(self):
        assert compute_case_name(12, 12, 0.0, 5) == "12x12J2=0D5"


class TestComputeWarmup:
    def test_small_lattice_clamps_to_minimum(self):
        assert compute_warmup(4, 4) == 100

    def test_medium_lattice_scales(self):
        assert compute_warmup(8, 8) == 128

    def test_large_lattice_clamps_to_maximum(self):
        assert compute_warmup(16, 16) == 500

    def test_rectangular(self):
        assert compute_warmup(12, 8) == 192


class TestComputeCgResidueRestart:
    def test_small_lattice(self):
        assert compute_cg_residue_restart(8, 8) == 10

    def test_large_lattice(self):
        assert compute_cg_residue_restart(16, 16) == 20

    def test_boundary(self):
        assert compute_cg_residue_restart(16, 8) == 20
        assert compute_cg_residue_restart(15, 15) == 10


def _default_params(**kwargs):
    defaults = dict(
        lx=8, ly=8, d=8, j2=0.0, bc="pbc",
        unit_lx=2, unit_ly=2,
        ntasks=56, optimizer="SR", lr=0.1,
        iterations=30,
        taus="0.5,0.2,0.1,0.05,0.02",
        step_cap=200,
        overrides={},
    )
    defaults.update(kwargs)
    return build_params(**defaults)


class TestBuildParams:
    def test_physics_unit_params(self):
        p = _default_params()
        pu = p["physics_unit"]
        assert pu["CaseParams"]["Lx"] == 2
        assert pu["CaseParams"]["Ly"] == 2
        assert pu["CaseParams"]["BoundaryCondition"] == "Periodic"

    def test_physics_full_params(self):
        p = _default_params()
        pf = p["physics_full"]
        assert pf["CaseParams"]["Lx"] == 8
        assert pf["CaseParams"]["Ly"] == 8

    def test_simple_update_has_tau_schedule(self):
        p = _default_params()
        su = p["simple_update"]["CaseParams"]
        assert su["TauScheduleEnabled"] is True
        assert su["TauScheduleTaus"] == "0.5,0.2,0.1,0.05,0.02"
        assert su["Dmin"] == 8
        assert su["Dmax"] == 8

    def test_vmc_trg_defaults(self):
        p = _default_params()
        vmc = p["vmc"]["CaseParams"]
        assert vmc["TRGDmin"] == 8
        assert vmc["TRGDmax"] == 16
        assert vmc["OptimizerType"] == "SR"

    def test_measure_doubles_mc_samples(self):
        p = _default_params()
        vmc_samples = p["vmc"]["CaseParams"]["MC_total_samples"]
        meas_samples = p["measure"]["CaseParams"]["MC_total_samples"]
        assert meas_samples == 2 * vmc_samples

    def test_override_applies(self):
        p = _default_params(overrides={"CGMaxIter": "500"})
        assert p["vmc"]["CaseParams"]["CGMaxIter"] == 500

    def test_obc_uses_bmps_not_trg(self):
        p = _default_params(bc="obc", unit_lx=4, unit_ly=4, lx=4, ly=4, d=4, ntasks=16)
        vmc = p["vmc"]["CaseParams"]
        assert "Dbmps_min" in vmc
        assert "TRGDmin" not in vmc


# ---- Slurm generation ----

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
            bc="pbc",
        )
        assert "simple_update" in content
        assert "sitps_tile" in content
        assert "--unit-ly 2" in content
        assert "--unit-lx 2" in content
        assert "source" in content
        assert "#SBATCH" in content

    def test_stage1_obc_omits_unit_flags(self):
        content = generate_slurm_stage1(
            job_name="hv_4x4D2_su",
            partition="256G56c",
            ntasks=1,
            walltime="00:15:00",
            case_dir="/share/home/wanghx/HeisenbergVMCPEPS/run/4x4J2=0D2/smoke_test",
            repo="/share/home/wanghx/HeisenbergVMCPEPS",
            unit_lx=4, unit_ly=4, full_lx=4, full_ly=4,
            bc="obc",
        )
        assert "simple_update" in content
        assert "sitps_tile" in content
        assert "--unit-ly" not in content
        assert "--unit-lx" not in content

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
        assert "cp" in content


# ---- File generation ----

class TestNewCaseFileGeneration:
    def test_new_case_creates_directory_structure(self, tmp_path):
        case_dir = tmp_path / "run" / "8x8J2=0D8" / "20260217_120000"
        params = _default_params()
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

        with open(case_dir / "params" / "vmc.json") as f:
            vmc = json.load(f)
        assert vmc["CaseParams"]["TRGDmin"] == 8


# ---- Error paths ----

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
                assert cmd_submit(args) == 1


class TestStatusErrorPaths:
    def test_status_ssh_failure(self):
        mock_result = MagicMock(returncode=1, stdout="", stderr="Connection timed out")
        with patch("cluster_run.ssh_run", return_value=mock_result):
            args = MagicMock()
            args.cluster = "susphy"
            args.user = "wanghx"
            assert cmd_status(args) == 1


# ---- Run journal ----

class TestParseEnergyTrajectory:
    def test_parses_csv(self, tmp_path):
        csv_file = tmp_path / "energy_trajectory.csv"
        csv_file.write_text(
            "iteration,energy,energy_err,gradient_norm\n"
            "0,-8.5,0.02,0.44\n"
            "1,-8.56,0.018,0.23\n"
            "2,-8.62,0.017,0.17\n"
        )
        result = parse_energy_trajectory(str(csv_file))
        assert result["final_energy"] == pytest.approx(-8.62)
        assert result["best_energy"] == pytest.approx(-8.62)
        assert result["best_iter"] == 2
        assert result["total_iterations"] == 3
        assert result["gradient_norm_final"] == pytest.approx(0.17)

    def test_missing_file_returns_none(self):
        assert parse_energy_trajectory("/nonexistent/path.csv") is None


class TestBuildJournalEntry:
    def test_builds_entry_from_meta_and_energy(self):
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
        assert result == 10000

    def test_ignores_different_config(self):
        journal = [
            {"params": {"lx": 4, "ly": 4, "D": 2, "bc": "obc"},
             "hpc": {"walltime_used_seconds": 100}},
        ]
        assert suggest_walltime(journal, 8, 8, 8, "pbc") is None


class TestFetchAppendsJournal:
    def test_fetch_creates_journal_entry(self, tmp_path):
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

        energy_dir = tmp_path / "remote_data" / "energy"
        energy_dir.mkdir(parents=True)
        (energy_dir / "energy_trajectory.csv").write_text(
            "iteration,energy,energy_err,gradient_norm\n"
            "0,-8.5,0.02,0.44\n"
            "1,-8.56,0.018,0.23\n"
            "2,-8.62,0.017,0.17\n"
        )

        journal_path = str(tmp_path / "run_journal.jsonl")

        def mock_scp(cluster, remote, local):
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
