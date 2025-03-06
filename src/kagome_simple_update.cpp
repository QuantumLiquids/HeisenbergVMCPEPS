// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-10-11
*
* Description: Simple Update for Heisenberg model in Kagome lattice.
*/


//#define PLAIN_TRANSPOSE 1

#include "./kagome_nn_on_sqr_peps_simple_update.h"
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
  Tensor ham_hei_tri = ham_hei_tri_terms[0] + ham_hei_tri_terms[1] + ham_hei_tri_terms[2];

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);
  // params.Ly is linear size of system in unit cell
  // peps size should be double of that.
  size_t peps_lx = 2 * params.Lx;
  size_t peps_ly = 2 * params.Ly;
  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, peps_ly, peps_lx);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    std::vector<std::vector<size_t>> activates(peps_ly, std::vector<size_t>(peps_lx, 0));
    size_t sz_int = 0;
    for (size_t y = 0; y < peps_ly - params.RemoveCorner; y++) {
      for (size_t x = 0; x < peps_lx - params.RemoveCorner; x++) {
        if ((x & 1) && (y & 1)) {
          activates[y][x] = 0;
        } else {
          activates[y][x] = sz_int % 2;
          sz_int += 1;
        }
      }
    }
    peps0.Initial(activates);
  }
  auto su_exe = new qlpeps::KagomeNNModelSquarePEPSSimpleUpdateExecutor<TenElemT, QNT>(update_para, peps0,
                                                                                       ham_hei_nn,
                                                                                       ham_hei_tri,
                                                                                       params.RemoveCorner);
  su_exe->Execute();
//  auto tps = qlpeps::TPS<TenElemT, QNT>(su_exe->GetPEPS());
//  tps.Dump();
  su_exe->DumpResult(peps_path, true);
  return 0;
}