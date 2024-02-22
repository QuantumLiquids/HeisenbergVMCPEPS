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

  Tensor ham_hei_tri_terms[3];
  for (size_t i = 0; i < 3; i++) {
    ham_hei_tri_terms[i] = Tensor({pb_in, pb_out, pb_in, pb_out, pb_in, pb_out});
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[0]({0, 0, 0, 0, i, i}) = 0.25;
    ham_hei_tri_terms[0]({1, 1, 1, 1, i, i}) = 0.25;
    ham_hei_tri_terms[0]({1, 1, 0, 0, i, i}) = -0.25;
    ham_hei_tri_terms[0]({0, 0, 1, 1, i, i}) = -0.25;
    ham_hei_tri_terms[0]({0, 1, 1, 0, i, i}) = 0.5;
    ham_hei_tri_terms[0]({1, 0, 0, 1, i, i}) = 0.5;
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[1]({0, 0, i, i, 0, 0}) = 0.25;
    ham_hei_tri_terms[1]({1, 1, i, i, 1, 1}) = 0.25;
    ham_hei_tri_terms[1]({1, 1, i, i, 0, 0}) = -0.25;
    ham_hei_tri_terms[1]({0, 0, i, i, 1, 1}) = -0.25;
    ham_hei_tri_terms[1]({0, 1, i, i, 1, 0}) = 0.5;
    ham_hei_tri_terms[1]({1, 0, i, i, 0, 1}) = 0.5;
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[2]({i, i, 0, 0, 0, 0}) = 0.25;
    ham_hei_tri_terms[2]({i, i, 1, 1, 1, 1}) = 0.25;
    ham_hei_tri_terms[2]({i, i, 1, 1, 0, 0}) = -0.25;
    ham_hei_tri_terms[2]({i, i, 0, 0, 1, 1}) = -0.25;
    ham_hei_tri_terms[2]({i, i, 0, 1, 1, 0}) = 0.5;
    ham_hei_tri_terms[2]({i, i, 1, 0, 0, 1}) = 0.5;
  }
  Tensor ham_hei_tri = 0.5 * ham_hei_tri_terms[0] + params.J2 * ham_hei_tri_terms[1] + 0.5 * ham_hei_tri_terms[2];

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);

//  size_t Ly = params.Ly, Lx = params.Lx;
  qlpeps::SquareLatticePEPS<TenElemT, U1QN> peps0(pb_out, params.Ly, params.Lx);
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
  if (params.J2 == 0) {
    auto su_exe = new qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, U1QN>(update_para, peps0,
                                                                                  ham_hei_nn);
    su_exe->Execute();
    auto tps = qlpeps::TPS<TenElemT, U1QN>(su_exe->GetPEPS());
    for (auto &tensor : tps) {
      tensor.Normalize();
    }
    tps.Dump();
    su_exe->DumpResult(peps_path, true);
  } else {
    auto su_exe = new qlpeps::SquareLatticeNNNSimpleUpdateExecutor<TenElemT, U1QN>(update_para, peps0,
                                                                                   ham_hei_nn,
                                                                                   ham_hei_tri);
    su_exe->Execute();
    auto tps = qlpeps::TPS<TenElemT, U1QN>(su_exe->GetPEPS());
    //TODO: if the first step vmc behave better, move into VMC package
    for (auto &tensor : tps) {
      tensor.Normalize();
    }
    tps.Dump();
    su_exe->DumpResult(peps_path, true);
  }

  return 0;
}