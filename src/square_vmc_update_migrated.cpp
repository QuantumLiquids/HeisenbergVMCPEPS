// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
* Migration Date: 2024-12-19
*
* Description: VMC Update for Heisenberg model - Migrated to VMCPEPSOptimizerExecutor.
*/

#include "./qldouble.h"
#include "params_parser.h"
#include "myutil.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer_params.h"
#include "qlpeps/optimizer/optimizer_params.h"

using namespace qlpeps;

using MCUpdater = MCUpdateSquareTNN3SiteExchange;

int main(int argc, char **argv) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  VMCUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;

  // Create separate parameter structures for the new API
  MonteCarloParams mc_params(
    params.MC_samples,
    // num_samples
    params.WarmUp,
    // num_warmup_sweeps
    params.MCLocalUpdateSweepsBetweenSample,
    // sweeps_between_samples
    "tps",
    // config_path
    Configuration(params.Ly, params.Lx) // alternative_init_config
  );

  // Use a default wavefunction path since VMCUpdateParams doesn't have it
  std::string wavefunction_path = "tps";

  PEPSParams peps_params(
    BMPSTruncatePara(params.Db_min,
                     params.Db_max,
                     params.TruncErr,
                     params.MPSCompressScheme,
                     std::make_optional<double>(params.TruncErr),
                     std::make_optional<size_t>(10)),
    wavefunction_path // wavefunction_path
  );

  OptimizerParams opt_params;
  // Fix step lengths assignment - use the first step length from the vector
  if (!params.step_len.empty()) {
    opt_params.core_params.step_lengths = {params.step_len[0], params.step_len[0], params.step_len[0]};
  } else {
    opt_params.core_params.step_lengths = {0.01, 0.01, 0.01}; // default values
  }
  opt_params.update_scheme = params.update_scheme;
  opt_params.cg_params = ConjugateGradientParams(params.CGMaxIter,
                                                 params.CGTol,
                                                 params.CGResidueRestart,
                                                 params.CGDiagShift);
  opt_params.core_params.max_iterations = params.step_len.size();
  VMCPEPSOptimizerParams optimize_para(opt_params, mc_params, peps_params);

  if (params.J2 == 0) {
    using Model = SquareSpinOneHalfXXZModel;
    VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    Model model; // Create default model instance

    if (IsFileExist(optimize_para.peps_params.wavefunction_path + "/tps_ten0_0_0.qlten")) {
      // Test if split index tps tensors exist
      executor = new VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                               params.Ly,
                                                                               params.Lx,
                                                                               comm,
                                                                               model);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                               tps,
                                                                               comm,
                                                                               model);
    }
    executor->Execute();
    delete executor;
  } else {
    using Model = SquareSpinOneHalfJ1J2XXZModel;
    VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.peps_params.wavefunction_path + "/tps_ten0_0_0.qlten")) {
      // Actually almost do the same thing
      executor = new VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                               params.Ly,
                                                                               params.Lx,
                                                                               comm,
                                                                               j1j2solver);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSOptimizerExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                               tps,
                                                                               comm,
                                                                               j1j2solver);
    }
    executor->Execute();
    delete executor;
  }
  MPI_Finalize();
  return 0;
}
