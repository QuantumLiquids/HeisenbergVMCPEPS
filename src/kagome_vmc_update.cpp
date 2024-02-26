//
// Created by haoxinwang on 12/10/2023.
//

#include "./qldouble.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps.h"
#include "kagome_hei_model_combined_tps_sample.h"
#include "spin_onehalf_heisenberg_kagome_model_sqrpeps_solver.h"
#include "./params_parser.h"
#include "myutil.h"

using namespace qlpeps;
using namespace std;

using TPSSampleT = KagomeCombinedTPSSampleLoaclFlip<TenElemT, U1QN>;

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

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

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

  qlpeps::VMCOptimizePara optimize_para;
  if (params.RemoveCorner) {
    Configuration init_config = GenerateInitialConfigurationInSmoothBoundary(params.Ly, params.Lx);
    optimize_para = qlpeps::VMCOptimizePara(
        BMPSTruncatePara(params.Db_min, params.Db_max,
                         params.TruncErr,
                         params.MPSCompressScheme,
                         std::make_optional<double>(params.TruncErr),
                         std::make_optional<size_t>(10)),
        params.MC_samples, params.WarmUp,
        params.MCLocalUpdateSweepsBetweenSample,
        init_config,
        params.step_len,
        params.update_scheme,
        ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));
  } else {
    optimize_para = qlpeps::VMCOptimizePara(
        BMPSTruncatePara(params.Db_min, params.Db_max,
                         params.TruncErr,
                         params.MPSCompressScheme,
                         std::make_optional<double>(params.TruncErr),
                         std::make_optional<size_t>(10)),
        params.MC_samples, params.WarmUp,
        params.MCLocalUpdateSweepsBetweenSample,
        occupation_num,
        params.Ly, params.Lx,
        params.step_len,
        params.update_scheme,
        ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));
  }

  using Model = KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver<TenElemT, U1QN>;
  VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
  Model kagome_heisenberg_model = Model(params.RemoveCorner);
  if (qlmps::IsPathExist(optimize_para.wavefunction_path)) {
    if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {//has split_index_tps
      executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                        params.Ly, params.Lx,
                                                                        world, kagome_heisenberg_model);
    } else {
      TPS<QLTEN_Double, U1QN> tps = TPS<QLTEN_Double, U1QN>(params.Ly, params.Lx);
      if (!tps.Load()) {
        std::cout << "Loading simple updated TPS files is broken." << std::endl;
        exit(-1);
      };
      executor = new VMCPEPSExecutor<QLTEN_Double, U1QN, TPSSampleT, Model>(optimize_para, tps,
                                                                            world, kagome_heisenberg_model);
    }
  } else {
    SquareLatticePEPS<QLTEN_Double, U1QN> peps(pb_out, 2 * params.Ly, 2 * params.Lx);
    if (!peps.Load(peps_path)) {
      std::cout << "Loading simple updated PEPS files is broken." << std::endl;
      exit(-2);
    };
    SplitIndexTPS<QLTEN_Double, U1QN> split_idx_tps = KagomeSquarePEPSToSplitIndexTPS(peps);
    if (!split_idx_tps.IsBondDimensionEven()) {
      std::cout << "Warning: Split Index TPS bond dimension  is not even!" << std::endl;
    }
    executor = new VMCPEPSExecutor<QLTEN_Double, U1QN, TPSSampleT, Model>(optimize_para, split_idx_tps,
                                                                          world, kagome_heisenberg_model);
  }

  executor->Execute();
  delete executor;

  return 0;
}
