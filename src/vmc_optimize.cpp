// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22 (simplified unified VMC entry)
*
* Description: Unified VMC Update with modern optimizer support for Heisenberg model.
*/

#include "./qldouble.h"
#include "enhanced_params_parser.h"
#include "qlpeps/qlpeps.h"
#include "model_updater_factory.h"
#include <fstream>
#include <variant>
#include "qlpeps/vmc_basic/spin_inversion_metadata.h"

using namespace qlpeps;

namespace {
const char *LBFGSStepModeToString(LBFGSStepMode mode) {
  switch (mode) {
    case LBFGSStepMode::kFixed: return "Fixed";
    case LBFGSStepMode::kStrongWolfe: return "StrongWolfe";
  }
  return "Unknown";
}

const char *MinSRSolverModeToString(MinSRSolverMode mode) {
  switch (mode) {
    case MinSRSolverMode::kAuto: return "Auto";
    case MinSRSolverMode::kReplicated: return "Replicated";
    case MinSRSolverMode::kDistributed: return "Distributed";
  }
  return "Unknown";
}
}  // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <physics_params.json> <vmc_algorithm_params.json>" << std::endl;
    std::cout << "Supported optimizers: SGD, Adam, AdaGrad, StochasticReconfiguration, LBFGS, MinSR" << std::endl;
    return -1;
  }
  
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  
  EnhancedVMCUpdateParams params(argv[1], argv[2]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.bmps_params.ThreadNum);

  if (rank == 0) {
    std::cout << "=== VMC Optimization ===" << std::endl;
    const auto &optimizer = params.optimizer_params;
    std::cout << "Optimizer: " << qlpeps::config::OptimizerTypeName(optimizer.algorithm_params)
              << std::endl;
    std::cout << "Learning Rate: " << optimizer.base_params.learning_rate << std::endl;
    std::cout << "Max Iterations: " << optimizer.base_params.max_iterations << std::endl;
    const auto &base = optimizer.base_params;
    if (base.lr_scheduler) {
      std::cout << "LR Scheduler: " << base.lr_scheduler->Name()
                << " (" << base.lr_scheduler->Describe() << ")" << std::endl;
    }
    if (base.clip_norm || base.clip_value) {
      std::cout << "Gradient Clipping: ";
      if (base.clip_norm) std::cout << "norm=" << *base.clip_norm << " ";
      if (base.clip_value) std::cout << "value=" << *base.clip_value << " ";
      std::cout << std::endl;
    }
    if (const auto *lbfgs = std::get_if<LBFGSParams>(&optimizer.algorithm_params)) {
      std::cout << "LBFGS Config: history=" << lbfgs->history_size
                << ", step_mode=" << LBFGSStepModeToString(lbfgs->step_mode)
                << ", max_eval=" << lbfgs->max_eval << std::endl;
      if (lbfgs->step_mode == LBFGSStepMode::kStrongWolfe) {
        std::cout << "  Strong-Wolfe: c1=" << lbfgs->wolfe_c1
                  << ", c2=" << lbfgs->wolfe_c2
                  << ", tol_grad=" << lbfgs->tolerance_grad
                  << ", tol_change=" << lbfgs->tolerance_change << std::endl;
      }
    }
    if (const auto *minsr = std::get_if<MinSRParams>(&optimizer.algorithm_params)) {
      std::cout << "MinSR Config: r_pinv=" << minsr->r_pinv
                << ", a_pinv=" << minsr->a_pinv
                << ", soft_cutoff=" << minsr->soft_cutoff
                << ", solver_mode=" << MinSRSolverModeToString(minsr->solver_mode)
                << std::endl;
    }
    const auto &initial = base.initial_step_selector;
    const auto &periodic = base.periodic_step_selector;
    if (initial.enabled || periodic.enabled) {
      std::cout << "Step Selectors:" << std::endl;
      std::cout << "  Initial: enabled=" << initial.enabled
                << ", max_line_search_steps=" << initial.max_line_search_steps
                << ", deterministic=" << initial.enable_in_deterministic << std::endl;
      std::cout << "  Periodic: enabled=" << periodic.enabled
                << ", every_n_steps=" << periodic.every_n_steps
                << ", phase_switch_ratio=" << periodic.phase_switch_ratio
                << ", deterministic=" << periodic.enable_in_deterministic << std::endl;
    }
    const auto &spike = optimizer.spike_recovery_params;
    if (spike.enable_auto_recover) {
      std::cout << "Spike Recovery: enabled (max_retries=" << spike.redo_mc_max_retries
                << ")" << std::endl;
    } else {
      std::cout << "Spike Recovery: disabled" << std::endl;
    }
    if (spike.enable_rollback) {
      std::cout << "  Rollback: enabled (sigma_k=" << spike.sigma_k << ")" << std::endl;
    }
    std::cout << "=================================" << std::endl;
  }

  // Initialize or load TPS/SITPS with unified basename
  std::string base = params.io_params.wavefunction_base; // default "tps"
  std::string tps_final = base + "final";
  // Note: we do not auto-fallback to lowest; user may manually copy lowest → final
  RequireSpinInversionMetadataCollectively(
      tps_final, MakeSpinInversionMetadata(params.spin_inversion_parity), comm);

  SplitIndexTPS<TenElemT, QNT> sitps;
  
  if (qlmps::IsPathExist(tps_final)) {
    if (rank == 0) std::cout << "Loading SplitIndexTPS from: " << tps_final << std::endl;
    // Debug-only probe: try load single-site tensor (0,0) first
    sitps = SplitIndexTPS<TenElemT, QNT>(params.physical_params.Ly, params.physical_params.Lx,
                                         params.physical_params.BoundaryCondition);
    sitps.Load(tps_final);
    if (sitps.GetBoundaryCondition() != params.physical_params.BoundaryCondition) {
      if (rank == 0) {
        std::cerr << "ERROR: BoundaryCondition mismatch between physics_params.json and loaded SplitIndexTPS.\n"
                  << "  physics BoundaryCondition = "
                  << ((params.physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic) ? "Periodic" : "Open")
                  << "\n  SITPS BoundaryCondition   = "
                  << ((sitps.GetBoundaryCondition() == qlpeps::BoundaryCondition::Periodic) ? "Periodic" : "Open")
                  << "\nPlease regenerate tpsfinal/ with the correct boundary condition." << std::endl;
      }
      MPI_Finalize();
      return -3;
    }
    if (rank == 0) std::cout << "Loaded SplitIndexTPS." << std::endl;
  } else {
    if (rank == 0) {
      std::cerr << "ERROR: Missing wavefunction directory '" << tps_final
                << "'. VMC requires SplitIndexTPS in tpsfinal/.\n"
                << "If you want to resume from the lowest snapshot, manually copy contents of tpslowest/ to tpsfinal/." << std::endl;
    }
    MPI_Finalize();
    return -2;
  }

  // Create and run VMC optimizer (backend-consistent dispatch).
  // TODO(MCRestrictU1): dispatch updater by params.mc_params.MCRestrictU1 as well.
  LogSamplerChoice(params.mc_params);
  if (rank == 0) {
    std::cout << "Starting optimization..." << std::endl;
  }
  RunVmcByModel<TenElemT, QNT>(params, sitps, comm, rank);
  if (rank == 0) {
    std::cout << "Optimization completed!" << std::endl;
  }

  MPI_Finalize();
  return 0;
}
