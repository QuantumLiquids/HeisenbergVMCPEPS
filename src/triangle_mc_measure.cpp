// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 15/01/2024.
*
* Description: Monte-Carlo measurement for Heisenberg model in triangle lattice.
*/

#include "gqpeps/algorithm/vmc_update/monte_carlo_measurement.h"
#include "gqpeps/algorithm/vmc_update/model_solvers/spin_onehalf_triangle_heisenberg_sqrpeps.h"
#include "gqpeps/algorithm/vmc_update/model_solvers/spin_onehalf_triangle_heisenbergJ1J2_sqrpeps.h"
#include "gqpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_3site_exchange.h"
#include "./gqdouble.h"
#include "./params_parser.h"
#include "myutil.h"

using namespace gqpeps;

using TPSSampleT = SquareTPSSample3SiteExchange<TenElemT, U1QN>;

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
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);

    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {// test if split index tps tensors exsit
      executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
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
    using Model = SpinOneHalfTriJ1J2HeisenbergSqrPEPS<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {// test if split index tps tensors exsit
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