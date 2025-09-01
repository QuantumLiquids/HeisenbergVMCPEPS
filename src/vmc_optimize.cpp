// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22 (simplified unified VMC entry)
*
* Description: Unified VMC Update with modern optimizer support for Heisenberg model.
*/

#include "./qldouble.h"
#include "enhanced_params_parser.h"
#include "myutil.h"
#include "qlpeps/qlpeps.h"
#include "qlpeps/api/vmc_api.h"
#include "model_updater_factory.h"
#include <fstream>

using namespace qlpeps;

using MCUpdater = MCUpdateSquareTNN3SiteExchange;
using Model = SquareSpinOneHalfXXZModel;

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <physics_params.json> <vmc_algorithm_params.json>" << std::endl;
    std::cout << "Supported optimizers: SGD, Adam, AdaGrad, StochasticReconfiguration" << std::endl;
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
    std::cout << "Optimizer: " << params.optimizer_type << std::endl;
    std::cout << "Learning Rate: " << params.learning_rate << std::endl;
    std::cout << "Max Iterations: " << params.max_iterations << std::endl;
    if (!params.lr_scheduler_type.empty()) {
      std::cout << "LR Scheduler: " << params.lr_scheduler_type << std::endl;
    }
    if (params.clip_norm || params.clip_value) {
      std::cout << "Gradient Clipping: ";
      if (params.clip_norm) std::cout << "norm=" << *params.clip_norm << " ";
      if (params.clip_value) std::cout << "value=" << *params.clip_value << " ";
      std::cout << std::endl;
    }
    std::cout << "=================================" << std::endl;
  }

  // Create VMC optimizer parameters (rank-aware config IO)
  auto vmc_params = params.CreateVMCOptimizerParams(rank);
  
  // Initialize or load TPS/SITPS with unified basename
  std::string base = params.wavefunction_base; // default "tps"
  std::string tps_final = base + "final";
  // Note: we do not auto-fallback to lowest; user may manually copy lowest → final

  SplitIndexTPS<TenElemT, QNT> sitps;
  
  if (qlmps::IsPathExist(tps_final)) {
    if (rank == 0) std::cout << "Loading SplitIndexTPS from: " << tps_final << std::endl;
    // Debug-only probe: try load single-site tensor (0,0) first
    sitps = SplitIndexTPS<TenElemT, QNT>(params.physical_params.Ly, params.physical_params.Lx);
    sitps.Load(tps_final);
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

  // Create and run VMC optimizer via high-level API wrapper
  // TODO(MCRestrictU1): dispatch MCUpdater by params.physical_params.MCRestrictU1
  LogSamplerChoice(params.physical_params);
  if (rank == 0) {
    std::cout << "Starting " << params.optimizer_type << " optimization..." << std::endl;
  }
  // Dispatch by model (SquareHeisenberg default; SquareHeisenbergJ1J2 fallback until wired)
  RunVmcByModel<TenElemT, QNT, MCUpdater>(params, sitps, comm);
  if (rank == 0) {
    std::cout << "Optimization completed!" << std::endl;
  }

  MPI_Finalize();
  return 0;
}


