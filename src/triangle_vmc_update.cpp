// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-28
*
* Description: VMC Update for Heisenberg model.
*/

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_triangle_heisenberg_sqrpeps.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_triangle_heisenbergJ1J2_sqrpeps.h"
#include "params_parser.h"
#include "myutil.h"

using namespace gqpeps;

using TPSSampleNNFlipT = SquareTPSSampleNNFlip<TenElemT, U1QN>;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  gqten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  gqpeps::VMCOptimizePara optimize_para(params.TruncErr, params.Db_min, params.Db_max,
                                        params.MPSCompressScheme,
                                        params.MC_samples, params.WarmUp,
                                        params.MCLocalUpdateSweepsBetweenSample,
                                        std::vector<size_t>{N / 2, N / 2},
                                        params.Ly, params.Lx,
                                        params.step_len,
                                        params.update_scheme);

  if (params.J2 == 0) {
    using Model = SpinOneHalfTriHeisenbergSqrPEPS<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT> *executor(nullptr);
    Model triangle_hei_solver;
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {
      executor = new VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT>(optimize_para,
                                                                              params.Ly, params.Lx,
                                                                              world, triangle_hei_solver);
    } else {
      TPS<TenElemT, U1QN> tps = TPS<TenElemT, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT>(optimize_para, tps,
                                                                              world, triangle_hei_solver);
    }
    executor->cg_params.max_iter = params.CGMaxIter;
    executor->cg_params.tolerance = params.CGTol;
    executor->cg_params.residue_restart_step = params.CGResidueRestart;
    executor->cg_params.diag_shift = params.CGDiagShift;
    executor->Execute();
    delete executor;
  } else {
    using Model = SpinOneHalfTriJ1J2HeisenbergSqrPEPS<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT> *executor(nullptr);
    double j2 = params.J2;
    Model trij1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) { //actually almostly do the same thing
      executor = new VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT>(optimize_para,
                                                                              params.Ly, params.Lx,
                                                                              world, trij1j2solver);
    } else {
      TPS<TenElemT, U1QN> tps = TPS<TenElemT, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, U1QN, Model, TPSSampleNNFlipT>(optimize_para, tps,
                                                                              world, trij1j2solver);
    }
    executor->cg_params.max_iter = params.CGMaxIter;
    executor->cg_params.tolerance = params.CGTol;
    executor->cg_params.residue_restart_step = params.CGResidueRestart;
    executor->cg_params.diag_shift = params.CGDiagShift;
    executor->Execute();
    delete executor;
  }

  return 0;
}
