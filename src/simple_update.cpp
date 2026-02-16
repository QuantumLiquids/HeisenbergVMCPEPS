// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: Unified Simple Update entry (Square Heisenberg/XY, Triangle Heisenberg).
*/

//#define PLAIN_TRANSPOSE 1

#include "qlpeps/algorithm/simple_update/square_lattice_nn_simple_update.h"
#include "qlpeps/algorithm/simple_update/square_lattice_nnn_simple_update.h"
#include "qlpeps/algorithm/simple_update/triangle_nn_on_sqr_peps_simple_update.h"
#include "qlpeps/qlpeps.h"
#include "qlpeps/api/conversions.h"
#include "./qldouble.h"
#include "./common_params.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

using StopReason = qlpeps::SimpleUpdateExecutor<TenElemT, QNT>::StopReason;

struct StageSummaryRecord {
  size_t index;
  double tau;
  size_t step_cap;
  bool converged;
  std::string stop_reason;
  size_t executed_steps;
};

std::string StopReasonToString(const StopReason stop_reason) {
  switch (stop_reason) {
    case StopReason::kMaxSteps:
      return "kMaxSteps";
    case StopReason::kAdvancedConverged:
      return "kAdvancedConverged";
    case StopReason::kNotRun:
    default:
      return "kNotRun";
  }
}

std::string FormatTauForPath(const double tau) {
  std::ostringstream oss;
  oss << std::defaultfloat << std::setprecision(15) << tau;
  std::string value = oss.str();
  const size_t exponent_pos = value.find_first_of("eE");
  if (exponent_pos == std::string::npos) {
    while (!value.empty() && value.back() == '0') {
      value.pop_back();
    }
    if (!value.empty() && value.back() == '.') {
      value.pop_back();
    }
  }
  if (value.empty()) {
    value = "0";
  }
  return value;
}

std::string BuildStageDirName(
    const size_t stage_index,
    const size_t stage_count,
    const double tau) {
  const size_t width = std::max<size_t>(2, std::to_string(stage_count).size());
  std::ostringstream oss;
  oss << "stage_" << std::setw(static_cast<int>(width)) << std::setfill('0') << stage_index
      << "_tau_" << FormatTauForPath(tau);
  return oss.str();
}

std::string EscapeJsonString(const std::string &input) {
  std::string escaped;
  escaped.reserve(input.size());
  for (const char c : input) {
    switch (c) {
      case '\\':
        escaped += "\\\\";
        break;
      case '"':
        escaped += "\\\"";
        break;
      case '\n':
        escaped += "\\n";
        break;
      default:
        escaped.push_back(c);
        break;
    }
  }
  return escaped;
}

void EnsureDirectory(const std::filesystem::path &dir_path) {
  std::error_code ec;
  std::filesystem::create_directories(dir_path, ec);
  if (ec) {
    throw std::runtime_error(
        "Failed to create directory '" + dir_path.string() + "': " + ec.message());
  }
}

void WriteScheduleSummaryJson(
    const std::string &dump_dir,
    const bool enabled,
    const size_t stage_count,
    const bool require_converged,
    const bool abort_on_stage_failure,
    const bool overall_success,
    const std::vector<StageSummaryRecord> &stages) {
  EnsureDirectory(dump_dir);
  const std::filesystem::path json_path =
      std::filesystem::path(dump_dir) / "schedule_summary.json";
  std::ofstream ofs(json_path);
  if (!ofs.is_open()) {
    throw std::runtime_error("Failed to open " + json_path.string() + " for writing.");
  }
  ofs << std::setprecision(15);
  ofs << "{\n"
      << "  \"enabled\": " << (enabled ? "true" : "false") << ",\n"
      << "  \"stage_count\": " << stage_count << ",\n"
      << "  \"require_converged\": " << (require_converged ? "true" : "false") << ",\n"
      << "  \"abort_on_stage_failure\": " << (abort_on_stage_failure ? "true" : "false") << ",\n"
      << "  \"overall_success\": " << (overall_success ? "true" : "false") << ",\n"
      << "  \"stages\": [\n";
  for (size_t i = 0; i < stages.size(); ++i) {
    const auto &stage = stages[i];
    ofs << "    {\n"
        << "      \"index\": " << stage.index << ",\n"
        << "      \"tau\": " << stage.tau << ",\n"
        << "      \"step_cap\": " << stage.step_cap << ",\n"
        << "      \"converged\": " << (stage.converged ? "true" : "false") << ",\n"
        << "      \"stop_reason\": \"" << EscapeJsonString(stage.stop_reason) << "\",\n"
        << "      \"executed_steps\": " << stage.executed_steps << "\n"
        << "    }";
    if (i + 1 < stages.size()) {
      ofs << ",";
    }
    ofs << "\n";
  }
  ofs << "  ]\n"
      << "}\n";
}

void WriteScheduleSummaryCsv(
    const std::string &dump_dir,
    const std::vector<StageSummaryRecord> &stages) {
  EnsureDirectory(dump_dir);
  const std::filesystem::path csv_path =
      std::filesystem::path(dump_dir) / "schedule_summary.csv";
  std::ofstream ofs(csv_path);
  if (!ofs.is_open()) {
    throw std::runtime_error("Failed to open " + csv_path.string() + " for writing.");
  }
  ofs << std::setprecision(15);
  ofs << "stage_index,tau,step_cap,converged,stop_reason,executed_steps\n";
  for (const auto &stage : stages) {
    ofs << stage.index << ","
        << stage.tau << ","
        << stage.step_cap << ","
        << (stage.converged ? "true" : "false") << ","
        << stage.stop_reason << ","
        << stage.executed_steps << "\n";
  }
}

void DumpSitpsAndPeps(
    qlpeps::SimpleUpdateExecutor<TenElemT, QNT> &su_exe,
    const std::string &sitps_dir,
    const std::string &peps_dir,
    const bool release_mem) {
  const auto &peps = su_exe.GetPEPS();
  auto tps = qlpeps::ToTPS<TenElemT, QNT>(peps);
  for (auto &tensor : tps) {
    tensor.Normalize();
  }
  qlpeps::SplitIndexTPS<TenElemT, QNT> sitps = qlpeps::ToSplitIndexTPS<TenElemT, QNT>(tps);
  sitps.Dump(sitps_dir);
  su_exe.DumpResult(peps_dir, release_mem);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <physics_params.json> <simple_update_algorithm_params.json>" << std::endl;
    return -1;
  }
  
  heisenberg_params::SimpleUpdateParams params(argv[1], argv[2]);

  // Model/data dispatch hints (Square XY handled by separate gate below)
  const std::string model = params.physical_params.ModelType;
  const bool is_triangle = (model.find("Triangle") != std::string::npos);
  const bool is_xy = (!is_triangle && (model == "SquareXY"));
  const auto bc = params.physical_params.BoundaryCondition;
  if (is_triangle && bc == qlpeps::BoundaryCondition::Periodic) {
    std::cerr << "ERROR: Triangle simple update with PBC is not supported." << std::endl;
    return -2;
  }

  Tensor did = Tensor({pb_in, pb_out});
  Tensor dsz = Tensor({pb_in, pb_out});
  Tensor dsp = Tensor({pb_in, pb_out});
  Tensor dsm = Tensor({pb_in, pb_out});
  Tensor ham_hei_nn = Tensor({pb_in, pb_out, pb_in, pb_out});

  did({0, 0}) = 1;
  did({1, 1}) = 1;
  dsz({0, 0}) = 0.5;
  dsz({1, 1}) = -0.5;
  dsp({0, 1}) = 1;
  dsm({1, 0}) = 1;

  // Heisenberg NN two-site gate
  ham_hei_nn({0, 0, 0, 0}) = 0.25;
  ham_hei_nn({1, 1, 1, 1}) = 0.25;
  ham_hei_nn({1, 1, 0, 0}) = -0.25;
  ham_hei_nn({0, 0, 1, 1}) = -0.25;
  ham_hei_nn({0, 1, 1, 0}) = 0.5;
  ham_hei_nn({1, 0, 0, 1}) = 0.5;

  // XY gate (SxSx + SySy) only has flip terms; no SzSz diagonal contribution
  Tensor ham_xy_nn = Tensor({pb_in, pb_out, pb_in, pb_out});
  ham_xy_nn({0, 1, 1, 0}) = 0.5; // |01><10|
  ham_xy_nn({1, 0, 0, 1}) = 0.5; // |10><01|

  // Select NN two-site Hamiltonian by ModelType
  Tensor ham_nn = is_xy ? ham_xy_nn : ham_hei_nn;
  auto ham_nnn = params.physical_params.J2 * ham_nn;

  // Triangle requires additional three-site term
  Tensor ham_hei_tri;
  if (is_triangle) {
    Tensor ham_hei_tri_terms[3];
    for (size_t i = 0; i < 3; i++) {
      ham_hei_tri_terms[i] = Tensor({pb_in, pb_out, pb_in, pb_out, pb_in, pb_out});
    }
    for (size_t i = 0; i < 2; i++) {
      ham_hei_tri_terms[0]({0, 0, 0, 0, i, i}) = 0.25;
      ham_hei_tri_terms[0]({1, 1, 1, 1, i, i}) = 0.25;
      ham_hei_tri_terms[0]({1, 1, 0, 0, i, i}) = -0.25;
      ham_hei_tri_terms[0]({0, 0, 1, 1, i, i}) = -0.25;
      ham_hei_tri_terms[0]({0, 1, 1, 0, i, i}) = 0.5;
      ham_hei_tri_terms[0]({1, 0, 0, 1, i, i}) = 0.5;
    }
    for (size_t i = 0; i < 2; i++) {
      ham_hei_tri_terms[1]({0, 0, i, i, 0, 0}) = 0.25;
      ham_hei_tri_terms[1]({1, 1, i, i, 1, 1}) = 0.25;
      ham_hei_tri_terms[1]({1, 1, i, i, 0, 0}) = -0.25;
      ham_hei_tri_terms[1]({0, 0, i, i, 1, 1}) = -0.25;
      ham_hei_tri_terms[1]({0, 1, i, i, 1, 0}) = 0.5;
      ham_hei_tri_terms[1]({1, 0, i, i, 0, 1}) = 0.5;
    }
    for (size_t i = 0; i < 2; i++) {
      ham_hei_tri_terms[2]({i, i, 0, 0, 0, 0}) = 0.25;
      ham_hei_tri_terms[2]({i, i, 1, 1, 1, 1}) = 0.25;
      ham_hei_tri_terms[2]({i, i, 1, 1, 0, 0}) = -0.25;
      ham_hei_tri_terms[2]({i, i, 0, 0, 1, 1}) = -0.25;
      ham_hei_tri_terms[2]({i, i, 0, 1, 1, 0}) = 0.5;
      ham_hei_tri_terms[2]({i, i, 1, 0, 0, 1}) = 0.5;
    }
    ham_hei_tri = ham_hei_tri_terms[0] + ham_hei_tri_terms[1] + ham_hei_tri_terms[2];
  }

  qlten::hp_numeric::SetTensorManipulationThreads(params.numerical_params.ThreadNum);

  const bool tau_schedule_enabled = params.tau_schedule.has_value();
  qlpeps::SimpleUpdatePara update_para = tau_schedule_enabled
      ? params.CreateSimpleUpdateParaForStage(
            params.tau_schedule.value().taus.front(),
            params.tau_schedule.value().step_caps.front())
      : params.CreateSimpleUpdatePara();

  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, params.physical_params.Ly, params.physical_params.Lx, bc);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    if (is_triangle) {
      // Three-sublattice ordered initial state for triangle lattice
      for (size_t y = 0; y < params.physical_params.Ly; y++) {
        for (size_t x = 0; x < params.physical_params.Lx; x++) {
          size_t sublattice_num;
          if (y >= x) {
            sublattice_num = (y - x) % 3;
          } else {
            sublattice_num = (3 - (x - y) % 3) % 3;
          }
          switch (sublattice_num) {
            case 0: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = 0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = 1; // spin down
                    break;
            case 1: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = -std::sqrt(3.0) / 2.0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
                    break;
            case 2: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = std::sqrt(3.0) / 2.0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
                    break;
          }
        }
      }
    } else {
      std::vector<std::vector<size_t>> activates(params.physical_params.Ly, std::vector<size_t>(params.physical_params.Lx));
      for (size_t y = 0; y < params.physical_params.Ly; y++) {
        for (size_t x = 0; x < params.physical_params.Lx; x++) {
          size_t sz_int = x + y;
          activates[y][x] = sz_int % 2;
        }
      }
      peps0.Initial(activates);
    }
  }

  std::unique_ptr<qlpeps::SimpleUpdateExecutor<TenElemT, QNT>> su_exe;
  if (is_triangle) {
    su_exe = std::make_unique<qlpeps::TriangleNNModelSquarePEPSSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_hei_nn, ham_hei_tri);
  } else if (std::abs(params.physical_params.J2) < 1e-15) {
    su_exe = std::make_unique<qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_nn);
  } else {
    su_exe = std::make_unique<qlpeps::SquareLatticeNNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_nn, ham_nnn);
  }

  const std::string sitps_final = "tpsfinal";
  int return_code = 0;
  if (!tau_schedule_enabled) {
    su_exe->Execute();
    if (params.advanced_stop.has_value()) {
      const auto &summary = su_exe->GetLastRunSummary();
      const std::string stop_reason = StopReasonToString(summary.stop_reason);
      std::cout << "Advanced stop summary: converged=" << std::boolalpha << summary.converged
                << ", stop_reason=" << stop_reason
                << ", executed_steps=" << summary.executed_steps
                << "/" << params.Step << std::endl;
    }
  } else {
    const auto &schedule = params.tau_schedule.value();
    std::vector<StageSummaryRecord> stage_summaries;
    stage_summaries.reserve(schedule.taus.size());
    bool overall_success = true;
    for (size_t stage = 0; stage < schedule.taus.size(); ++stage) {
      const size_t stage_index = stage + 1;
      const double stage_tau = schedule.taus[stage];
      const size_t stage_step_cap = schedule.step_caps[stage];
      std::cout << "=== Tau stage " << stage_index << "/" << schedule.taus.size()
                << ": tau=" << stage_tau
                << ", step_cap=" << stage_step_cap
                << " ===" << std::endl;

      su_exe->update_para = params.CreateSimpleUpdateParaForStage(stage_tau, stage_step_cap);
      su_exe->Execute();
      const auto &summary = su_exe->GetLastRunSummary();
      const std::string stop_reason = StopReasonToString(summary.stop_reason);
      std::cout << "=== Stage result: converged=" << std::boolalpha << summary.converged
                << ", executed_steps=" << summary.executed_steps
                << ", stop_reason=" << stop_reason
                << " ===" << std::endl;

      stage_summaries.push_back(StageSummaryRecord{
          stage_index,
          stage_tau,
          stage_step_cap,
          summary.converged,
          stop_reason,
          summary.executed_steps
      });

      if (schedule.dump_each_stage) {
        const std::filesystem::path stage_dir =
            std::filesystem::path(schedule.dump_dir) /
            BuildStageDirName(stage_index, schedule.taus.size(), stage_tau);
        EnsureDirectory(stage_dir);
        DumpSitpsAndPeps(
            *su_exe,
            (stage_dir / "tpsfinal").string(),
            (stage_dir / "peps").string(),
            false);
      }

      const bool stage_failed = schedule.require_converged && !summary.converged;
      if (stage_failed) {
        overall_success = false;
        if (schedule.abort_on_stage_failure) {
          return_code = -3;
          break;
        }
      }
    }

    WriteScheduleSummaryJson(
        schedule.dump_dir,
        true,
        schedule.taus.size(),
        schedule.require_converged,
        schedule.abort_on_stage_failure,
        overall_success,
        stage_summaries);
    WriteScheduleSummaryCsv(schedule.dump_dir, stage_summaries);
  }

  DumpSitpsAndPeps(*su_exe, sitps_final, peps_path, true);

  std::cout << "Simple Update completed." << std::endl;
  std::cout << "SplitIndexTPS saved to: " << sitps_final << std::endl;
  std::cout << "PEPS saved to: " << peps_path << std::endl;

  return return_code;
}
