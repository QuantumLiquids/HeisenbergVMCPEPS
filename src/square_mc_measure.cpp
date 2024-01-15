//
// Created by haoxinwang on 15/01/2024.
//


#include "gqpeps/algorithm/vmc_update/monte_carlo_measurement.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_heisenberg_square.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_squareJ1J2.h"
#include "./gqdouble.h"
#include "./params_parser.h"
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
    using Model = SpinOneHalfHeisenbergSquare<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleNNFlipT, Model> *executor(nullptr);

    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleNNFlipT, Model>(optimize_para,
                                                                                            params.Ly, params.Lx,
                                                                                            world);
    } else {

    }

    if (params.ReplicaTest) {
      executor->ReplicaTest();
    } else {
      executor->Execute();
    }
    delete executor;
    std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);
  } else {
    using Model = SpinOneHalfJ1J2HeisenbergSquare<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleNNFlipT, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleNNFlipT, Model>(optimize_para,
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