// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: Simple Update for Heisenberg model.
*/

//#define PLAIN_TRANSPOSE 1

#include "qlpeps/algorithm/simple_update/square_lattice_nn_simple_update.h"
#include "qlpeps/algorithm/simple_update/square_lattice_nnn_simple_update.h"
#include "./qldouble.h"
#include "./params_parser.h"

int main(int argc, char **argv) {
  SimpleUpdateParams params(argv[1]);

  Tensor did = Tensor({pb_in, pb_out});
  Tensor dsz = Tensor({pb_in, pb_out});
  Tensor dsp = Tensor({pb_in, pb_out});
  Tensor dsm = Tensor({pb_in, pb_out});
  Tensor ham_hei_nn = Tensor({pb_in, pb_out, pb_in, pb_out});

  did({0, 0}) = 1;
  did({1, 1}) = 1;
  dsz({0, 0}) = 0.5;
  dsz({1, 1}) = -0.5;
  dsp({0, 1}) = 1;
  dsm({1, 0}) = 1;

  ham_hei_nn({0, 0, 0, 0}) = 0.25;
  ham_hei_nn({1, 1, 1, 1}) = 0.25;
  ham_hei_nn({1, 1, 0, 0}) = -0.25;
  ham_hei_nn({0, 0, 1, 1}) = -0.25;
  ham_hei_nn({0, 1, 1, 0}) = 0.5;
  ham_hei_nn({1, 0, 0, 1}) = 0.5;

  auto ham_hei_nnn = params.J2 * ham_hei_nn;

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);

//  size_t Ly = params.Ly, Lx = params.Lx;
  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, params.Ly, params.Lx);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    std::vector<std::vector<size_t>> activates(params.Ly, std::vector<size_t>(params.Lx));
    for (size_t y = 0; y < params.Ly; y++) {
      for (size_t x = 0; x < params.Lx; x++) {
        size_t sz_int = x + y;
        activates[y][x] = sz_int % 2;
      }
    }
    peps0.Initial(activates);
  }
  std::unique_ptr<qlpeps::SimpleUpdateExecutor<TenElemT, QNT>> su_exe;
  if (params.J2 == 0) {
      su_exe = std::make_unique<qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_hei_nn);
  } else {
      su_exe = std::make_unique<qlpeps::SquareLatticeNNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_hei_nn, ham_hei_nnn);
  }
  su_exe->Execute();
  auto tps = qlpeps::TPS<TenElemT, QNT>(su_exe->GetPEPS());
  for (auto &tensor : tps) {
      tensor.Normalize();
  }
  tps.Dump();
  su_exe->DumpResult(peps_path, true);

  return 0;
}