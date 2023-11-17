//
// Created by haoxinwang on 04/11/2023.
//

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "./params_parser.h"
#include "kagome_mc_measure.h"

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
                                        params.MCLocalUpdateSweepsBetweenSample,
                                        occupation_num, params.step_len,
                                        params.update_scheme);
  optimize_para.mc_sweep_scheme = CompressedLatticeKagomeLocalUpdate;


  using Model = KagomeSpinOneHalfHeisenbergMeasurementSolver<TenElemT, U1QN>;
  KagomeMeasurementExecutor<TenElemT, U1QN, Model> *executor(nullptr);
  executor = new KagomeMeasurementExecutor<TenElemT, U1QN, Model>(optimize_para,
                                                                  params.Ly, params.Lx,
                                                                  world);
  bool replica_test = params.Continue_from_VMC;
  if (replica_test) {
    executor->ReplicaTest();
  } else {
    executor->Execute();
  }

//  executor->cg_params.max_iter = params.CGMaxIter;
//  executor->cg_params.tolerance = params.CGTol;
//  executor->cg_params.residue_restart_step = params.CGResidueRestart;
//  executor->cg_params.diag_shift = params.CGDiagShift;

  delete executor;
  return 0;
}
