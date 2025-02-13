// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: VMC Update for Heisenberg model.
*/

#include "./qldouble.h"
#include "params_parser.h"
#include "myutil.h"
#include "qlpeps/qlpeps.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSample3SiteExchange<TenElemT, QNT>;

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
    using Model = SpinOneHalfHeisenbergSquare<TenElemT, QNT>;
    VMCPEPSExecutor<TenElemT, QNT, TPSSampleT, Model> *executor(nullptr);

    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exist
      executor = new VMCPEPSExecutor<TenElemT, QNT, TPSSampleT, Model>(optimize_para,
                                                                        params.Ly, params.Lx,
                                                                        comm);
    } else {
      TPS<QLTEN_Double, QNT> tps = TPS<QLTEN_Double, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<TenElemT, QNT, TPSSampleT, Model>(optimize_para, tps,
                                                                        comm);
    }
    executor->Execute();
    delete executor;
  } else {
    using Model = SpinOneHalfJ1J2HeisenbergSquare<QLTEN_Double, QNT>;
    VMCPEPSExecutor<QLTEN_Double, QNT, TPSSampleT, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) { //actually almost do the same thing
      executor = new VMCPEPSExecutor<QLTEN_Double, QNT, TPSSampleT, Model>(optimize_para,
                                                                            params.Ly, params.Lx,
                                                                            comm, j1j2solver);
    } else {
      TPS<QLTEN_Double, QNT> tps = TPS<QLTEN_Double, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<QLTEN_Double, QNT, TPSSampleT, Model>(optimize_para, tps,
                                                                            comm, j1j2solver);
    }
    executor->Execute();
    delete executor;
  }
  return 0;
}
