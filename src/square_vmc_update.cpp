// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: VMC Update for Heisenberg model.
*/

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_heisenberg_square.h"    // SpinOneHalfHeisenbergSquare
#include "gqpeps/algorithm/vmc_update/model_energy_solvers/spin_onehalf_squareJ1J2.h"           // SpinOneHalfJ1J2HeisenbergSquare
#include "params_parser.h"
#include "myutil.h"

using namespace gqpeps;

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
    VMCPEPSExecutor<TenElemT, U1QN, Model> *executor(nullptr);

    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) {// test if split index tps tensors exsit
      executor = new VMCPEPSExecutor<TenElemT, U1QN, Model>(optimize_para,
                                                            params.Ly, params.Lx,
                                                            world);
    } else {
      TPS<GQTEN_Double, U1QN> tps = TPS<GQTEN_Double, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para, tps,
                                                                world);
    }

    executor->Execute();
    delete executor;
  } else {
    using Model = SpinOneHalfJ1J2HeisenbergSquare<GQTEN_Double, U1QN>;
    VMCPEPSExecutor<GQTEN_Double, U1QN, Model> *executor(nullptr);
    double j2 = params.J2;
    Model j1j2solver(j2);
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.gqten")) { //actually almostly do the same thing
      executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para,
                                                                params.Ly, params.Lx,
                                                                world, j1j2solver);
    } else {
      TPS<GQTEN_Double, U1QN> tps = TPS<GQTEN_Double, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-2);
      };
      executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para, tps,
                                                                world, j1j2solver);
    }

    executor->Execute();
    delete executor;
  }
  return 0;
}
