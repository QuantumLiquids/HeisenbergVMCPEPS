//
// Created by haoxinwang on 04/11/2023.
//

#include "./gqdouble.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "./params_parser.h"
#include "kagome_mc_measure.h"

using namespace gqpeps;

Configuration GenerateInitialConfigurationInSmoothBoundary(size_t ly, size_t lx) {
  size_t sys_ly = 2 * ly, sys_lx = 2 * lx;
  std::vector<std::vector<size_t>> activates(sys_ly, std::vector<size_t>(sys_lx, 0));
  size_t sz_int = 0;
  for (size_t y = 0; y < sys_ly - 1; y++) {
    for (size_t x = 0; x < sys_lx - 1; x++) {
      if ((x & 1) && (y & 1)) {
        activates[y][x] = 0;
      } else {
        activates[y][x] = sz_int % 2;
        sz_int += 1;
      }
    }
  }
  Configuration config(ly, lx);
  for (size_t row = 0; row < ly; row++) {
    for (size_t col = 0; col < lx; col++) {
      config({row, col}) = activates[2 * row][2 * col]
                           + 2 * activates[2 * row + 1][2 * col]
                           + 4 * activates[2 * row][2 * col + 1];
    }
  }
  return config;
}


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

  gqpeps::VMCOptimizePara optimize_para;
  if (params.RemoveCorner) {
    Configuration init_config = GenerateInitialConfigurationInSmoothBoundary(params.Ly, params.Lx);
    optimize_para = gqpeps::VMCOptimizePara(params.TruncErr, params.Db_min, params.Db_max,
                                            params.MPSCompressScheme,
                                            params.MC_samples, params.WarmUp,
                                            params.MCLocalUpdateSweepsBetweenSample,
                                            init_config,
                                            params.step_len,
                                            params.update_scheme);
  } else {
    optimize_para = gqpeps::VMCOptimizePara(params.TruncErr, params.Db_min, params.Db_max,
                                            params.MPSCompressScheme,
                                            params.MC_samples, params.WarmUp,
                                            params.MCLocalUpdateSweepsBetweenSample,
                                            occupation_num,
                                            params.Ly, params.Lx,
                                            params.step_len,
                                            params.update_scheme);
  }
  optimize_para.mc_sweep_scheme = CompressedLatticeKagomeLocalUpdate;


  using Model = KagomeSpinOneHalfHeisenbergMeasurementSolver<TenElemT, U1QN>;
  KagomeMeasurementExecutor<TenElemT, U1QN, Model> *executor(nullptr);
  executor = new KagomeMeasurementExecutor<TenElemT, U1QN, Model>(optimize_para,
                                                                  params.Ly, params.Lx,
                                                                  world);
  if (params.ReplicaTest) {
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
