#!/usr/bin/env python3
"""Small driver-level MPI regression; all state and output files live in a temporary directory."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

# Set by each repository's copy.
MODEL = "heisenberg"


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    args = parser.parse_args()
    binary = args.binary.resolve()
    tj = MODEL == "tj"
    physics = (dict(Lx=3, Ly=2, Hole=2, t=1.0, t2=0.0, J=0.3, V=0.1, Mu=0.2)
               if tj else dict(Lx=3, Ly=2, J2=0.0, RemoveCorner=False,
                               ModelType="SquareHeisenberg", BoundaryCondition="Open"))
    algorithm = dict(Dmin=1, Dmax=2, TruncErr=1e-10, Tau=0.05, Step=1, ThreadNum=1,
                     AdvancedStopEnabled=False)
    with tempfile.TemporaryDirectory(prefix=f"{MODEL}-su-mpi-") as temporary:
        root = Path(temporary)

        def run(name, ranks, phys, algo, initial=None, pinning=False, failure=False):
            work = root / name
            work.mkdir()
            for filename, params in (("physics.json", phys), ("algorithm.json", algo)):
                (work / filename).write_text(json.dumps({"CaseParams": params}))
            if initial is not None:
                shutil.copytree(initial / "peps", work / "peps")
            command = [str(binary), "physics.json", "algorithm.json"]
            if pinning:
                (work / "pinning.json").write_text(json.dumps({"CaseParams": {
                    "X": [0, 2], "Y": [0, 1], "Add_Mu": [0.1, -0.05], "AFM_Pinning": 0.02}}))
                command.append("pinning.json")
            if ranks > 1 or name == "launched_singleton":
                command = [args.mpiexec, "-n", str(ranks)] + command
            result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=40)
            output = result.stdout + result.stderr
            check((result.returncode != 0) == failure,
                  f"{name}: return code {result.returncode}\n{output}")
            if not failure:
                check(output.count("Simple Update completed.") == 1, f"{name}: duplicate completion output")
                check(any((work / "peps").iterdir()), f"{name}: no PEPS output")
                check(any((work / "tpsfinal").iterdir()), f"{name}: no TPS output")
                check(("SimpleUpdate MPI process grid:" in output) == (ranks > 1),
                      f"{name}: wrong serial/MPI dispatch")
            return work, output

        # Preserve legacy serial execution and create one authoritative (possibly random) input.
        seed, _ = run("seed", 1, physics, algorithm)
        run("launched_singleton", 1, physics, algorithm, seed)
        schedule = dict(algorithm, Step=3, TauScheduleEnabled=True,
                        TauScheduleTaus="0.05,0.02", TauScheduleStepCaps="3,3",
                        TauScheduleRequireConverged=True, TauScheduleDumpEachStage=True,
                        TauScheduleDumpDir="stages", TauScheduleAbortOnStageFailure=True,
                        AdvancedStopEnabled=True, AdvancedStopEnergyAbsTol=1e6,
                        AdvancedStopEnergyRelTol=1e6, AdvancedStopLambdaRelTol=1e6,
                        AdvancedStopPatience=1, AdvancedStopMinSteps=2)
        for nnn in (False, True):
            phys = dict(physics)
            phys["t2" if tj else "J2"] = 0.2 if nnn else 0.0
            # Also cover a periodic spin system and the t-J pinning constructor.
            if nnn and not tj:
                phys["BoundaryCondition"] = "Periodic"
                periodic_seed, _ = run("periodic_seed", 1, phys, algorithm)
                initial = periodic_seed
            else:
                initial = seed
            results = []
            for ranks in (2, 4):
                work, output = run(f"{'nnn' if nnn else 'nn'}_{ranks}", ranks, phys,
                                   schedule, initial, pinning=tj and nnn)
                summary = json.loads((work / "stages/schedule_summary.json").read_text())
                check(summary["overall_success"] and len(summary["stages"]) == 2,
                      "Schedule did not complete both stages")
                check(all(s["converged"] and s["executed_steps"] == 2 for s in summary["stages"]),
                      "Collective advanced stop disagrees with the expected step count")
                check(len(list((work / "stages").glob("stage_*/peps"))) == 2,
                      "Missing stage checkpoints")
                results.append(work)
            # Same input and layered order must yield the same persisted tensors across rank counts.
            for directory in ("peps", "tpsfinal"):
                left = {p.relative_to(results[0] / directory): p.read_bytes()
                        for p in (results[0] / directory).rglob("*") if p.is_file()}
                right = {p.relative_to(results[1] / directory): p.read_bytes()
                         for p in (results[1] / directory).rglob("*") if p.is_file()}
                check(left == right, f"{directory}: persisted state differs between 2 and 4 ranks")
            run(f"restart_{nnn}", 2, phys, algorithm, results[1], pinning=tj and nnn)

        # Preserve each driver's existing nonconverged-stage policy and nonzero exit status.
        failed = dict(schedule, TauScheduleStepCaps="1,1", AdvancedStopMinSteps=10)
        work, _ = run("stage_failure", 2, physics, failed, seed, failure=True)
        summary = json.loads((work / "stages/schedule_summary.json").read_text())
        check(not summary["overall_success"], "Failed stage reported success")
        check(len(summary["stages"]) == (2 if tj else 1), "Stage-failure policy changed")
        if not tj:
            _, output = run("triangle_rejected", 2, dict(physics, ModelType="TriangleHeisenberg"),
                            algorithm, failure=True)
            check("square NN/NNN" in output, "Missing unsupported-model diagnostic")
    print(f"PASS: {MODEL} serial/MPI dispatch, NN/NNN, 2/4 ranks, checkpoints, restart and stopping")


if __name__ == "__main__":
    main()
