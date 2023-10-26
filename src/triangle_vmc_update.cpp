// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-28
*
* Description: VMC Update for Heisenberg model.
*/

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_triangle_heisenberg_sqrpeps.h"
#include "params_parser.h"


using namespace gqpeps;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  gqten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  gqpeps::VMCOptimizePara optimize_para(params.TruncErr, params.Db_min, params.Db_max,
                                        params.MC_samples, params.WarmUp,
                                        params.MCLocalUpdateSweepsBetweenSample,
                                        {N / 2, N / 2}, params.step_len,
                                        params.update_scheme);


  using Model = SpinOneHalfTriHeisenbergSqrPEPS<GQTEN_Double, U1QN>;
  VMCPEPSExecutor<GQTEN_Double, U1QN, Model> *executor(nullptr);
  Model triangle_hei_solver;
  if (params.Continue_from_VMC) {
    executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para,
                                                              params.Ly, params.Lx,
                                                              world, triangle_hei_solver);
  } else {
    TPS<GQTEN_Double, U1QN> tps = TPS<GQTEN_Double, U1QN>(params.Ly, params.Lx);
    if (!tps.Load()) {
      std::cout << "Loading simple updated TPS files is broken." << std::endl;
      exit(-2);
    };
    executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para, tps,
                                                              world, triangle_hei_solver);
  }

  executor->Execute();
  delete executor;

  return 0;
}
