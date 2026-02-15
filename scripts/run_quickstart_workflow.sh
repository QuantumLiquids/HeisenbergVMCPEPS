#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run full PEPS workflow and keep reproducible artifacts for later analysis.

Workflow:
  simple_update -> vmc_optimize -> mc_measure -> plot

Usage:
  scripts/run_quickstart_workflow.sh [options]

Options:
  --profile <name>            Preset: local22_pbc (default) | cluster44_obc
  --build-dir <path>          Build directory containing binaries (default: build_llvm_sdk)
  --run-root <path>           Root directory for new run folders (default: run)
  --tag <text>                Extra suffix in run directory name
  --physics <path>            Physics params JSON (override preset)
  --simple-update <path>      Simple update params JSON (override preset)
  --vmc <path>                VMC params JSON (override preset)
  --measure <path>            Measurement params JSON (override preset)
  --mpi-ranks <int>           MPI ranks for vmc/measure (default from profile)
  --oversubscribe             Add --oversubscribe to mpirun (useful on local testing)
  --plot-title <text>         Title for plot/workflow/plot_energy_trajectory.py
  --no-plot                   Skip plotting step
  --help                      Show this message

Examples:
  scripts/run_quickstart_workflow.sh
  scripts/run_quickstart_workflow.sh --profile cluster44_obc --mpi-ranks 16
  scripts/run_quickstart_workflow.sh --profile cluster44_obc --mpi-ranks 16 --oversubscribe
  scripts/run_quickstart_workflow.sh --physics /path/to/physics.json --vmc /path/to/vmc.json
EOF
}

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

abs_path() {
  local path="$1"
  local base="$2"
  if [[ "$path" = /* ]]; then
    echo "$path"
  else
    echo "$base/$path"
  fi
}

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "[error] Missing file: $path" >&2
    exit 2
  fi
}

require_executable() {
  local path="$1"
  if [[ ! -x "$path" ]]; then
    echo "[error] Missing executable or not executable: $path" >&2
    exit 2
  fi
}

quote_cmd() {
  printf "%q " "$@"
}

run_step() {
  local step_name="$1"
  shift
  local log_file="$run_dir/${step_name}.log"

  echo "[$(timestamp)] STEP ${step_name}: $(quote_cmd "$@")" >> "$run_dir/commands.log"

  set +e
  (
    cd "$run_dir"
    "$@"
  ) >"$log_file" 2>&1
  local rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "FAILED step=${step_name} rc=${rc}" > "$run_dir/status.txt"
    echo "[error] Step '${step_name}' failed with exit code ${rc}. See: $log_file" >&2
    echo "[error] Last 120 lines:" >&2
    tail -n 120 "$log_file" >&2 || true
    exit "$rc"
  fi
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

profile="local22_pbc"
build_dir="build_llvm_sdk"
run_root="run"
tag=""
physics=""
simple_update=""
vmc=""
measure=""
mpi_ranks=""
oversubscribe=0
do_plot=1
plot_title=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      profile="${2:-}"
      shift 2
      ;;
    --build-dir)
      build_dir="${2:-}"
      shift 2
      ;;
    --run-root)
      run_root="${2:-}"
      shift 2
      ;;
    --tag)
      tag="${2:-}"
      shift 2
      ;;
    --physics)
      physics="${2:-}"
      shift 2
      ;;
    --simple-update)
      simple_update="${2:-}"
      shift 2
      ;;
    --vmc)
      vmc="${2:-}"
      shift 2
      ;;
    --measure)
      measure="${2:-}"
      shift 2
      ;;
    --mpi-ranks)
      mpi_ranks="${2:-}"
      shift 2
      ;;
    --oversubscribe)
      oversubscribe=1
      shift
      ;;
    --plot-title)
      plot_title="${2:-}"
      shift 2
      ;;
    --no-plot)
      do_plot=0
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "[error] Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

case "$profile" in
  local22_pbc)
    default_physics="params/quickstart/physics_local_2x2_pbc.json"
    default_su="params/quickstart/simple_update_local_2x2_pbc.json"
    default_vmc="params/quickstart/vmc_local_2x2_pbc_n1.json"
    default_measure="params/quickstart/measure_local_2x2_pbc_n1.json"
    default_mpi=1
    default_title="Square Heisenberg 2x2 PBC (quick local)"
    ;;
  cluster44_obc)
    default_physics="params/quickstart/physics_cluster_4x4_obc.json"
    default_su="params/quickstart/simple_update_cluster_4x4_obc.json"
    default_vmc="params/quickstart/vmc_cluster_4x4_obc_n16.json"
    default_measure="params/quickstart/measure_cluster_4x4_obc_n16.json"
    default_mpi=16
    default_title="Square Heisenberg 4x4 OBC (cluster)"
    ;;
  *)
    echo "[error] Unsupported profile: $profile" >&2
    echo "Supported profiles: local22_pbc, cluster44_obc" >&2
    exit 2
    ;;
esac

if [[ -z "$physics" ]]; then physics="$default_physics"; fi
if [[ -z "$simple_update" ]]; then simple_update="$default_su"; fi
if [[ -z "$vmc" ]]; then vmc="$default_vmc"; fi
if [[ -z "$measure" ]]; then measure="$default_measure"; fi
if [[ -z "$mpi_ranks" ]]; then mpi_ranks="$default_mpi"; fi
if [[ -z "$plot_title" ]]; then plot_title="$default_title"; fi

build_dir="$(abs_path "$build_dir" "$repo_root")"
run_root="$(abs_path "$run_root" "$repo_root")"
physics="$(abs_path "$physics" "$repo_root")"
simple_update="$(abs_path "$simple_update" "$repo_root")"
vmc="$(abs_path "$vmc" "$repo_root")"
measure="$(abs_path "$measure" "$repo_root")"

require_file "$physics"
require_file "$simple_update"
require_file "$vmc"
require_file "$measure"

simple_update_bin="$build_dir/simple_update"
vmc_bin="$build_dir/vmc_optimize"
measure_bin="$build_dir/mc_measure"
plot_script="$repo_root/plot/workflow/plot_energy_trajectory.py"

require_executable "$simple_update_bin"
require_executable "$vmc_bin"
require_executable "$measure_bin"

if ! [[ "$mpi_ranks" =~ ^[0-9]+$ ]] || [[ "$mpi_ranks" -lt 1 ]]; then
  echo "[error] --mpi-ranks must be a positive integer, got: $mpi_ranks" >&2
  exit 2
fi

run_stamp="$(date +%Y%m%d_%H%M%S)"
run_name="${run_stamp}_${profile}"
if [[ -n "$tag" ]]; then
  run_name="${run_name}_${tag}"
fi
run_dir="$run_root/$run_name"

mkdir -p "$run_dir/params_used"

cp "$physics" "$run_dir/params_used/physics.json"
cp "$simple_update" "$run_dir/params_used/simple_update.json"
cp "$vmc" "$run_dir/params_used/vmc.json"
cp "$measure" "$run_dir/params_used/measure.json"

if command -v shasum >/dev/null 2>&1; then
  (
    cd "$run_dir/params_used"
    shasum -a 256 physics.json simple_update.json vmc.json measure.json > checksums.sha256
  )
fi

{
  echo "repo_root=$repo_root"
  echo "run_dir=$run_dir"
  echo "profile=$profile"
  echo "build_dir=$build_dir"
  echo "physics=$physics"
  echo "simple_update=$simple_update"
  echo "vmc=$vmc"
  echo "measure=$measure"
  echo "mpi_ranks=$mpi_ranks"
  echo "oversubscribe=$oversubscribe"
  echo "plot_enabled=$do_plot"
  echo "plot_title=$plot_title"
  echo "start_time=$(timestamp)"
  echo "host=$(hostname)"
  echo "user=$(whoami)"
  if command -v git >/dev/null 2>&1; then
    echo "git_head=$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || true)"
    echo "git_branch=$(git -C "$repo_root" rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
  fi
} > "$run_dir/run_info.txt"

echo "[info] run_dir: $run_dir"
echo "[info] profile: $profile"
echo "[info] mpi_ranks: $mpi_ranks"

run_step "01_simple_update" "$simple_update_bin" "$physics" "$simple_update"

if [[ "$mpi_ranks" -eq 1 ]]; then
  run_step "02_vmc" "$vmc_bin" "$physics" "$vmc"
  run_step "03_measure" "$measure_bin" "$physics" "$measure"
else
  mpirun_cmd=(mpirun)
  if [[ "$oversubscribe" -eq 1 ]]; then
    mpirun_cmd+=(--oversubscribe)
  fi
  mpirun_cmd+=(-n "$mpi_ranks")
  run_step "02_vmc" "${mpirun_cmd[@]}" "$vmc_bin" "$physics" "$vmc"
  run_step "03_measure" "${mpirun_cmd[@]}" "$measure_bin" "$physics" "$measure"
fi

if [[ "$do_plot" -eq 1 ]]; then
  if [[ ! -f "$plot_script" ]]; then
    echo "[warn] Plot script not found: $plot_script" | tee -a "$run_dir/warnings.log"
  elif [[ ! -f "$run_dir/energy/energy_trajectory.csv" ]]; then
    echo "[warn] Missing energy CSV: $run_dir/energy/energy_trajectory.csv" | tee -a "$run_dir/warnings.log"
  else
    run_step "04_plot" \
      python3 "$plot_script" \
      --csv "$run_dir/energy/energy_trajectory.csv" \
      --out "$run_dir/energy/energy_trajectory.png" \
      --title "$plot_title"
  fi
fi

find "$run_dir" -type f | sed "s|^$run_dir/||" | sort > "$run_dir/files_manifest.txt"
find "$run_dir" -type f | sed "s|^$run_dir/||" | \
  grep -E '(^|/)stats/.+\.csv$|(^|/)samples/psi\.csv$|one_point_functions\.csv$|two_point_functions\.csv$|(^|/)energy_sample_data/' \
  > "$run_dir/measurement_outputs.txt" || true

{
  echo "=== CHECKPOINTS ==="
  grep -nE "Simple Update completed\\.|SplitIndexTPS saved to: tpsfinal|Loading SplitIndexTPS from: tpsfinal|Configuration validation|Optimization completed!|MONTE-CARLO MEASUREMENT PROGRAM|Rank 0: statistic data finished\\." \
    "$run_dir/01_simple_update.log" "$run_dir/02_vmc.log" "$run_dir/03_measure.log" || true
  echo
  echo "=== ARTIFACTS ==="
  ls -ld "$run_dir"/tpsfinal "$run_dir"/peps "$run_dir"/tpslowest "$run_dir"/energy \
    "$run_dir"/stats "$run_dir"/samples 2>/dev/null || true
  ls -l "$run_dir"/energy/energy_trajectory.csv "$run_dir"/energy/energy_trajectory.png 2>/dev/null || true
} > "$run_dir/checkpoints_and_artifacts.txt"

if [[ -f "$run_dir/energy/energy_trajectory.csv" ]]; then
  head -n 20 "$run_dir/energy/energy_trajectory.csv" > "$run_dir/energy_head20.csv"
fi

echo "SUCCESS end_time=$(timestamp)" > "$run_dir/status.txt"
echo "[ok] Workflow finished."
echo "[ok] Run directory: $run_dir"
echo "[ok] Key files:"
echo "  - $run_dir/run_info.txt"
echo "  - $run_dir/checkpoints_and_artifacts.txt"
echo "  - $run_dir/measurement_outputs.txt"
echo "  - $run_dir/files_manifest.txt"
