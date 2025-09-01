// SPDX-License-Identifier: LGPL-3.0-only

//
// Monte Carlo measurement entry using the new PEPS API (qlpeps/api/vmc_api.h)
// - No legacy MCMeasurementPara / Executor types
// - Model dispatch by J2
// - Prefer SplitIndexTPS from tpsfinal/; fallback to splitting TPS
//

#include "qlpeps/qlpeps.h"
#include "qlpeps/api/vmc_api.h"
#include "./qldouble.h"
#include "enhanced_measure_params_parser.h"
#include "model_updater_factory.h"
#include "myutil.h"

using namespace qlpeps;

using MCUpdater = MCUpdateSquareTNN3SiteExchange;

int main(int argc, char **argv) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  if (argc != 3) {
    if (rank == 0) std::cout << "Usage: " << argv[0] << " <physics_params.json> <measure_algorithm_params.json>" << std::endl;
    MPI_Finalize();
    return -1;
  }

  EnhancedMCMeasureParams params(argv[1], argv[2]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.bmps_params.ThreadNum);

  // Build MC measurement params (new API)
  const size_t N = params.physical_params.Lx * params.physical_params.Ly;
  Configuration init_config(params.physical_params.Ly, params.physical_params.Lx);
  bool warmed_up = false;
  {
    // Try to load warmed-up configuration (rank-aware) from tpsfinal/
    std::string load_dir = "tpsfinal";
    if (qlmps::IsPathExist(load_dir)) {
      warmed_up = init_config.Load(load_dir, (size_t) rank);
    }
    if (!warmed_up) {
      // Random initialization with U1 occupancy: half up, half down
      OccupancyNum occupancy(2, 0);
      occupancy[1] = N / 2;                 // spin up
      occupancy[0] = N - occupancy[1];      // spin down
      init_config.Random(occupancy);
    }
  }

  MonteCarloParams mc_params_obj(
      params.mc_params.MC_samples,
      params.mc_params.WarmUp,
      params.mc_params.MCLocalUpdateSweepsBetweenSample,
      init_config,
      warmed_up,              // warmed-up if loaded
      "tpsfinal"             // dump final configuration for reuse
  );
  PEPSParams peps_params_obj(params.CreateBMPSPara());
  MCMeasurementParams measurement_params(mc_params_obj, peps_params_obj, "./");

  // TODO(MCRestrictU1): dispatch MCUpdater by params.physical_params.MCRestrictU1
  LogSamplerChoice(params.physical_params);

  // Load SplitIndexTPS from tpsfinal/ if exists; otherwise split from TPS
  SplitIndexTPS<TenElemT, QNT> sitps;
  const std::string tps_final_dir = "tpsfinal";
  if (qlmps::IsPathExist(tps_final_dir)) {
    if (rank == 0) std::cout << "Loading SplitIndexTPS from: " << tps_final_dir << std::endl;
    sitps = SplitIndexTPS<TenElemT, QNT>(params.physical_params.Ly, params.physical_params.Lx);
    sitps.Load(tps_final_dir);
  } else {
    if (rank == 0) std::cout << "SplitIndexTPS not found. Loading TPS and splitting indices..." << std::endl;
    TPS<QLTEN_Double, QNT> tps(params.physical_params.Ly, params.physical_params.Lx);
    if (!tps.Load()) {
      if (rank == 0) std::cerr << "ERROR: Failed to load TPS from current directory." << std::endl;
      MPI_Finalize();
      return -2;
    }
    sitps = SplitIndexTPS<TenElemT, QNT>(tps);
  }

  // Dispatch by model (same policy as VMC)
  RunMeasureByModel<TenElemT, QNT, MCUpdater>(params.physical_params, measurement_params, sitps, comm);

  MPI_Finalize();
  return 0;
}

