// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-02-10
*
* Description: VMC Update for Heisenberg model in triangle lattice, without conservation on Sz quantum number.
*/

#include "./qldouble.h"
#include "qlpeps/qlpeps.h"
#include "params_parser.h"
#include "myutil.h"

using namespace qlpeps;

using MCUpdater = MCUpdateSquareNNFullSpaceUpdate;

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

  if (params.J2 == 0) {
    using Model = SpinOneHalfTriHeisenbergSqrPEPS<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    Model triangle_hei_solver;
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {
      executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                      params.Ly, params.Lx,
                                                                      world, triangle_hei_solver);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para, tps,
                                                                      world, triangle_hei_solver);
    }
    executor->Execute();
    delete executor;
  } else {
    using Model = SpinOneHalfTriJ1J2HeisenbergSqrPEPS<TenElemT, QNT>;
    VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    double j2 = params.J2;
    Model trij1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) { //actually almostly do the same thing
      executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para,
                                                                      params.Ly, params.Lx,
                                                                      world, trij1j2solver);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, QNT, MCUpdater, Model>(optimize_para, tps,
                                                                      world, trij1j2solver);
    }
    executor->Execute();
    delete executor;
  }

  return 0;
}
