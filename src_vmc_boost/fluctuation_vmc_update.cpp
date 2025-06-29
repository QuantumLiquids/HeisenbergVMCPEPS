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

  using Model = SpinOneHalfWaveFunctionFluctuationSqrPEPS<TenElemT, QNT>;
  VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
  Model triangle_hei_solver;
  if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {
    executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                    params.Ly, params.Lx,
                                                                    comm, triangle_hei_solver);
  } else {
    TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
    if (!tps.Load()) {
      std::cout << "Loading simple updated TPS files is broken." << std::endl;
      exit(-2);
    };
    executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para, tps,
                                                                    comm, triangle_hei_solver);
  }
  executor->Execute();
  delete executor;
  return 0;
}
