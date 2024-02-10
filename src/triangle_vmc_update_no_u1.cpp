// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-02-10
*
* Description: VMC Update for Heisenberg model in triangle lattice, without conservation on Sz quantum number.
*/

#include "./qldouble.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_triangle_heisenberg_sqrpeps.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_triangle_heisenbergJ1J2_sqrpeps.h"
#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_full_space_nn_flip.h"
#include "params_parser.h"
#include "myutil.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSampleFullSpaceNNFlip<TenElemT, U1QN>;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  qlpeps::VMCOptimizePara optimize_para(
      BMPSTruncatePara(params.Db_min, params.Db_max,
                       params.TruncErr,
                       params.MPSCompressScheme),
      params.MC_samples, params.WarmUp,
      params.MCLocalUpdateSweepsBetweenSample,
      std::vector<size_t>{N / 2, N / 2},
      params.Ly, params.Lx,
      params.step_len,
      params.update_scheme,
      ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));

  if (params.J2 == 0) {
    using Model = SpinOneHalfTriHeisenbergSqrPEPS<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
    Model triangle_hei_solver;
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {
      executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                        params.Ly, params.Lx,
                                                                        world, triangle_hei_solver);
    } else {
      TPS<TenElemT, U1QN> tps = TPS<TenElemT, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para, tps,
                                                                        world, triangle_hei_solver);
    }
    executor->Execute();
    delete executor;
  } else {
    using Model = SpinOneHalfTriJ1J2HeisenbergSqrPEPS<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
    double j2 = params.J2;
    Model trij1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) { //actually almostly do the same thing
      executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                        params.Ly, params.Lx,
                                                                        world, trij1j2solver);
    } else {
      TPS<TenElemT, U1QN> tps = TPS<TenElemT, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para, tps,
                                                                        world, trij1j2solver);
    }
    executor->Execute();
    delete executor;
  }

  return 0;
}
