//
// Created by haoxinwang on 15/01/2024.
//


#include "qlpeps/algorithm/vmc_update/monte_carlo_measurement.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_heisenberg_square.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_squareJ1J2.h"
#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_3site_exchange.h"
#include "./qldouble.h"
#include "./params_parser.h"
#include "myutil.h"

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
      params.update_scheme);

  if (params.J2 == 0) {
    using Model = SpinOneHalfHeisenbergSquare<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);

    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                                      params.Ly, params.Lx,
                                                                                      world);
    } else {
      TPS<QLTEN_Double, U1QN> tps = TPS<QLTEN_Double, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                                      SplitIndexTPS<TenElemT,
                                                                                                    U1QN>(tps),
                                                                                      world);
    }

    if (params.ReplicaTest) {
      executor->ReplicaTest();
    } else {
      executor->Execute();
    }
    delete executor;
    std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);
  } else {
    using Model = SpinOneHalfPlanarJ1J2HeisenbergSquare<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                                      params.Ly, params.Lx,
                                                                                      world, j1j2solver);
    } else {

    }

    if (params.ReplicaTest) {
      executor->ReplicaTest();
    } else {
      executor->Execute();
    }
  }

  return 0;
}