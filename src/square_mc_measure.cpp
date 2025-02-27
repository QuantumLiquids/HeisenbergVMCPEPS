//
// Created by haoxinwang on 15/01/2024.
//


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
    using Model = SquareSpinOneHalfXXZModel;
    MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);

    if (IsFileExist(
        measurement_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    params.Ly, params.Lx,
                                                                                    comm);
    } else {
      TPS<QLTEN_Double, QNT> tps = TPS<QLTEN_Double, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    SplitIndexTPS<TenElemT,
                                                                                                  QNT>(tps),
                                                                                    comm);
    }

    executor->Execute();
    executor->OutputEnergy();
    delete executor;
    std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);
  } else {
    using Model = SquareSpinOneHalfJ1J2XXZModel;
    MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(
        measurement_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    params.Ly, params.Lx,
                                                                                    comm, j1j2solver);
    } else {
      TPS<QLTEN_Double, QNT> tps = TPS<QLTEN_Double, QNT>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new MonteCarloMeasurementExecutor<TenElemT, QNT, MCUpdater, Model>(measurement_para,
                                                                                    SplitIndexTPS<TenElemT,
                                                                                                  QNT>(tps),
                                                                                    comm, j1j2solver);
    }
    executor->Execute();
    executor->OutputEnergy();
  }
  MPI_Finalize();
  return 0;
}