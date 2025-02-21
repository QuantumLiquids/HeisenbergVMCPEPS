// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 15/01/2024.
*
* Description: Monte-Carlo measurement for Heisenberg model in triangle lattice.
*/

#include "qlpeps/qlpeps.h"
#include "./qldouble.h"
#include "./params_parser.h"
#include "myutil.h"

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
  qlpeps::MCMeasurementPara measurement_para(
      BMPSTruncatePara(params.Db_min, params.Db_max,
                       params.TruncErr,
                       params.MPSCompressScheme,
                       std::make_optional<double>(params.TruncErr),
                       std::make_optional<size_t>(10)),
      params.MC_samples, params.WarmUp,
      params.MCLocalUpdateSweepsBetweenSample,
      std::vector<size_t>{N / 2, N / 2},
      params.Ly, params.Lx);

  if (params.J2 == 0) {
    using Model = SpinOneHalfTriHeisenbergSqrPEPS;
    MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);

    if (IsFileExist(
        measurement_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    params.Ly, params.Lx,
                                                                                    comm);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      SplitIndexTPS<TenElemT, QNT> split_index_tps(tps);
      split_index_tps.NormalizeAllSite();
      std::cout << "Show Split Index TPS tensors Info :" << std::endl;
      for (auto &tens : split_index_tps) {
        for (auto &ten : tens) {
          ten.ConciseShow();
        }
      }
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    split_index_tps,
                                                                                    comm);
    }

    executor->Execute();
    executor->OutputEnergy();
    delete executor;
    std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);
  } else {
    using Model = SpinOneHalfTriJ1J2HeisenbergSqrPEPS;
    MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(
        measurement_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    params.Ly, params.Lx,
                                                                                    comm, j1j2solver);
    } else {
      TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      SplitIndexTPS<TenElemT, QNT> split_index_tps(tps);
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    split_index_tps,
                                                                                    comm, j1j2solver);
    }

    executor->Execute();
    executor->OutputEnergy();
  }
  MPI_Finalize();
  return 0;
}