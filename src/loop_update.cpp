// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2026-02-20
*
* Description: Square Heisenberg loop-update driver with simple-update-like parameter surface.
*/

#include "qlpeps/algorithm/loop_update/loop_update.h"
#include "qlpeps/api/conversions.h"
#include "qlpeps/qlpeps.h"
#include "./qldouble.h"
#include "./common_params.h"
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace {

using LoopGateT = qlpeps::LoopGates<Tensor>;
using StopReason = qlpeps::LoopUpdateExecutor<TenElemT, QNT>::StopReason;

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

#ifdef U1SYM
IndexT vb_out = IndexT({QNSctT(QNT(0), 5)}, qlten::TenIndexDirType::OUT);
#else
IndexT vb_out = IndexT({QNSctT(qn0, 5)}, qlten::TenIndexDirType::OUT);
#endif
IndexT vb_in = InverseIndex(vb_out);

LoopGateT GenerateSquareHeisenbergLoopGates(
    const double tau,
    const size_t n0,
    const size_t n1,
    const size_t n2,
    const size_t n3) {
  const std::array<size_t, 4> ns = {n0, n1, n2, n3};
  LoopGateT gates;
  for (size_t i = 0; i < 4; ++i) {
    gates[i] = Tensor({vb_in, pb_in, pb_out, vb_out});
    Tensor &gate = gates[i];
    // A->A (0->0): identity
    gate({0, 0, 0, 0}) = 1.0;
    gate({0, 1, 1, 0}) = 1.0;
    // A->Cz (0->1): Sz
    gate({0, 0, 0, 1}) = 0.5;
    gate({0, 1, 1, 1}) = -0.5;
    // A->C+ (0->2): S+
    gate({0, 1, 0, 2}) = 1.0;
    // A->C- (0->3): S-
    gate({0, 0, 1, 3}) = 1.0;
    // Cz->B (1->4): -tau/n * Sz
    gate({1, 0, 0, 4}) = -0.5 * tau / static_cast<double>(ns[i]);
    gate({1, 1, 1, 4}) = 0.5 * tau / static_cast<double>(ns[i]);
    // C+->B (2->4): -tau/(2n) * S-
    gate({2, 0, 1, 4}) = -tau / (2.0 * static_cast<double>(ns[i]));
    // C-->B (3->4): -tau/(2n) * S+
    gate({3, 1, 0, 4}) = -tau / (2.0 * static_cast<double>(ns[i]));
    // B->A (4->0): identity closure
    gate({4, 0, 0, 0}) = 1.0;
    gate({4, 1, 1, 0}) = 1.0;
  }
  return gates;
}

void GenerateSquareHeisenbergAllEvolveGates(
    const double tau,
    const size_t ly,
    const size_t lx,
    const qlpeps::BoundaryCondition bc,
    qlpeps::DuoMatrix<LoopGateT> &evolve_gates) {
  if (bc == qlpeps::BoundaryCondition::Periodic) {
    for (size_t row = 0; row < ly; ++row) {
      for (size_t col = 0; col < lx; ++col) {
        evolve_gates({row, col}) = GenerateSquareHeisenbergLoopGates(tau, 2, 2, 2, 2);
      }
    }
    return;
  }

  for (size_t row = 0; row < ly - 1; ++row) {
    for (size_t col = 0; col < lx - 1; ++col) {
      const size_t n_left = (col == 0) ? 1 : 2;
      const size_t n_upper = (row == 0) ? 1 : 2;
      const size_t n_right = (col == lx - 2) ? 1 : 2;
      const size_t n_lower = (row == ly - 2) ? 1 : 2;
      evolve_gates({row, col}) =
          GenerateSquareHeisenbergLoopGates(tau, n_left, n_upper, n_right, n_lower);
    }
  }
}

void DumpSitpsAndPeps(
    qlpeps::LoopUpdateExecutor<TenElemT, QNT> &loop_exe,
    const std::string &sitps_dir,
    const std::string &peps_dir,
    const bool release_mem) {
  const auto &peps = loop_exe.GetPEPS();
  auto tps = qlpeps::ToTPS<TenElemT, QNT>(peps);
  for (auto &tensor : tps) {
    tensor.Normalize();
  }
  qlpeps::SplitIndexTPS<TenElemT, QNT> sitps = qlpeps::ToSplitIndexTPS<TenElemT, QNT>(tps);
  sitps.Dump(sitps_dir);
  loop_exe.DumpResult(peps_dir, release_mem);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " <physics_params.json> <loop_update_algorithm_params.json>" << std::endl;
    return -1;
  }

  try {
    heisenberg_params::LoopUpdateParams params(argv[1], argv[2]);
    const auto bc = params.physical_params.BoundaryCondition;
    const size_t ly = params.physical_params.Ly;
    const size_t lx = params.physical_params.Lx;

    if (bc == qlpeps::BoundaryCondition::Open && (ly < 2 || lx < 2)) {
      throw std::invalid_argument(
          "loop_update on OBC requires Lx>=2 and Ly>=2 so at least one plaquette exists.");
    }

    qlten::hp_numeric::SetTensorManipulationThreads(params.numerical_params.ThreadNum);

    const size_t gate_rows = (bc == qlpeps::BoundaryCondition::Periodic) ? ly : ly - 1;
    const size_t gate_cols = (bc == qlpeps::BoundaryCondition::Periodic) ? lx : lx - 1;
    qlpeps::DuoMatrix<LoopGateT> evolve_gates(gate_rows, gate_cols);
    GenerateSquareHeisenbergAllEvolveGates(params.Tau, ly, lx, bc, evolve_gates);

    if (!qlmps::IsPathExist(peps_path)) {
      throw std::runtime_error(
          "loop_update requires an existing PEPS input at '" + peps_path +
          "'. Run simple_update first (or provide peps/ in current working directory).");
    }

    qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, ly, lx, bc);
    if (!peps0.Load(peps_path)) {
      throw std::runtime_error("Failed to load PEPS from: " + peps_path);
    }

    const auto loop_para = params.CreateLoopUpdatePara();
    auto loop_exe = std::make_unique<qlpeps::LoopUpdateExecutor<TenElemT, QNT>>(
        loop_para, evolve_gates, peps0);
    loop_exe->Execute();

    if (params.advanced_stop.has_value()) {
      const auto &summary = loop_exe->GetLastRunSummary();
      const std::string stop_reason = StopReasonToString(summary.stop_reason);
      std::cout << "Advanced stop summary: converged=" << std::boolalpha << summary.converged
                << ", stop_reason=" << stop_reason
                << ", executed_steps=" << summary.executed_steps
                << "/" << params.Step << std::endl;
    }

    const std::string sitps_final = "tpsfinal";
    DumpSitpsAndPeps(*loop_exe, sitps_final, peps_path, true);

    std::cout << "Loop Update completed." << std::endl;
    std::cout << "SplitIndexTPS saved to: " << sitps_final << std::endl;
    std::cout << "PEPS saved to: " << peps_path << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return -2;
  }
}
