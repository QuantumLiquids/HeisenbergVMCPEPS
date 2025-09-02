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

// Local helper: init/load half-up-half-down Configuration with U1 occupancy
static inline std::pair<Configuration, bool> InitOrLoadConfigHalfU1(
    const heisenberg_params::PhysicalParams &physical_params,
    const std::string &configuration_load_dir,
    int rank) {
  const size_t N = physical_params.Lx * physical_params.Ly;
  const size_t spin_up_sites = N / 2;
  OccupancyNum occupancy(2, 0);
  occupancy[1] = spin_up_sites;
  occupancy[0] = N - spin_up_sites;
  Configuration config(physical_params.Ly, physical_params.Lx, occupancy);
  bool warmed_up = false;
  if (!configuration_load_dir.empty()) {
    warmed_up = config.Load(configuration_load_dir, static_cast<size_t>(rank));
  }
  return {config, warmed_up};
}

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
  auto [init_config, warmed_up] = InitOrLoadConfigHalfU1(
      params.physical_params,
      std::string("tpsfinal"),
      rank);

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

