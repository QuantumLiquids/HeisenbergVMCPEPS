// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-11-01
*
* Description: Minimize the error of estimation on the wave function component
*/

#include "../src/qldouble.h"
#include "qlpeps/qlpeps.h"
#include "../src/params_parser.h"
#include "../src/myutil.h"
#include "./spin_onehalf_wave_function_fluctuation_sqrpeps.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSample3SiteExchange<TenElemT, U1QN>;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  qlpeps::VMCOptimizePara optimize_para(
      BMPSTruncatePara(params.Db_min, params.Db_max,
                       params.TruncErr,
                       params.MPSCompressScheme,
                       std::make_optional<double>(params.TruncErr),
                       std::make_optional<size_t>(10)),
      params.MC_samples, params.WarmUp,
      params.MCLocalUpdateSweepsBetweenSample,
      std::vector<size_t>{N / 2, N / 2},
      params.Ly, params.Lx,
      params.step_len,
      params.update_scheme,
      ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));

  using Model = SpinOneHalfWaveFunctionFluctuationSqrPEPS<TenElemT, U1QN>;
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
  return 0;
}
