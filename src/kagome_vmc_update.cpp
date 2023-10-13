//
// Created by haoxinwang on 12/10/2023.
//

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "spin_onehalf_heisenberg_kagome_model_sqrpeps_solver.h"
#include "./params_parser.h"

using namespace gqpeps;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  gqten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  std::vector<size_t> occupation_num(8, 0);

  size_t occupation = N / 6;  //no 0 and 7;
  size_t remain_occupation = N - 6 * occupation;
  for (size_t i = 1; i < 7; i++) {
    occupation_num[i] = occupation;
  }
  if (remain_occupation > 0) {
    occupation_num[1] += remain_occupation / 2;
    occupation_num[3] += remain_occupation / 2;
  }
  if (remain_occupation % 2 == 1) {
    std::cout << "lattice site number = odd, not support." << std::endl;
    exit(-1);
  }

  gqpeps::VMCOptimizePara optimize_para(params.TruncErr, params.Db_min, params.Db_max,
                                        params.MC_samples, params.WarmUp,
                                        occupation_num, params.step_len,
                                        params.update_scheme);
  optimize_para.mc_sweep_sheme = CompressedLatticeKagomeLocalUpdate;


  using Model = KagomeSpinOneHalfHeisenbergSquare<TenElemT, U1QN>;
  VMCPEPSExecutor<TenElemT, U1QN, Model> *executor(nullptr);
  if (params.Continue_from_VMC) {
    executor = new VMCPEPSExecutor<TenElemT, U1QN, Model>(optimize_para,
                                                          params.Ly, params.Lx,
                                                          world);
  } else {
    SquareLatticePEPS<GQTEN_Double, U1QN> peps(pb_out, 2 * params.Ly, 2 * params.Lx);
    if (!peps.Load(peps_path)) {
      std::cout << "Loading simple updated PEPS files is broken." << std::endl;
      exit(-2);
    };
    SplitIndexTPS<GQTEN_Double, U1QN> split_idx_tps = KagomeSquarePEPSToSplitIndexTPS(peps);

    executor = new VMCPEPSExecutor<GQTEN_Double, U1QN, Model>(optimize_para, split_idx_tps,
                                                              world);
  }

  executor->Execute();
  delete executor;

  return 0;
}
