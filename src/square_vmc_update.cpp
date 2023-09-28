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
#include "gqmps2/case_params_parser.h"

struct VMCUpdateParams : public gqmps2::CaseParamsParserBasic {
  VMCUpdateParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    J2 = ParseDouble("J2");
    Db_min = ParseInt("Dbmps_min");
    Db_max = ParseInt("Dbmps_max");
    TruncErr = ParseDouble("TruncErr");
    MC_samples = ParseInt("MC_samples");
    WarmUp = ParseInt("WarmUp");
    Continue_from_VMC = ParseBool("Continue_from_VMC");
    size_t update_times = ParseInt("UpdateNum");
    step_len = std::vector<double>(update_times);
    if (update_times > 0) {
      step_len[0] = ParseDouble("StepLengthFirst");
      double step_len_change = ParseDouble("StepLengthDecrease");
      for (size_t i = 1; i < update_times; i++) {
        step_len[i] = step_len[0] - i * step_len_change;
      }
    }
    update_scheme = (gqpeps::WAVEFUNCTION_UPDATE_SCHEME) ParseInt("UpdateScheme");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Ly;
  size_t Lx;
  double J2;
  size_t Db_min;
  size_t Db_max;
  size_t TruncErr;
  size_t MC_samples;
  size_t WarmUp;
  bool Continue_from_VMC;
  gqpeps::WAVEFUNCTION_UPDATE_SCHEME update_scheme;
  std::vector<double> step_len;
  size_t ThreadNum;
};

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
                                        {N / 2, N / 2}, params.step_len,
                                        params.update_scheme);

  if (params.J2 == 0) {
    using Model = SpinOneHalfHeisenbergSquare<TenElemT, U1QN>;
    VMCPEPSExecutor<TenElemT, U1QN, Model> *executor(nullptr);

    if (params.Continue_from_VMC) {
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
    if (params.Continue_from_VMC) {
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
