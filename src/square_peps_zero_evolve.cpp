// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: Simple Update with tau = 0
*/

//#define PLAIN_TRANSPOSE 1

#include "qlpeps/algorithm/simple_update/square_lattice_nn_simple_update.h"
#include "./qldouble.h"
#include "./params_parser.h"

using namespace qlpeps;
int main(int argc, char **argv) {
  SimpleUpdateParams params(argv[1]);
  Tensor ham_hei_nn = Tensor({pb_in, pb_out, pb_in, pb_out});
  Tensor ham_hei_tri = Tensor({pb_in, pb_out, pb_in, pb_out, pb_in, pb_out});

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, 1,
                                       params.Dmax, params.Dmax,
                                       params.TruncErr);

  SplitIndexTPS<TenElemT, QNT> split_index_tps(params.Ly, params.Lx);
  if (!split_index_tps.Load()) {
    exit(1);
  }
  TPS<TenElemT, QNT> tps = split_index_tps.GroupIndices(pb_out);

//  size_t Ly = params.Ly, Lx = params.Lx;
  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, params.Ly, params.Lx);
  peps0.Gamma = tps;
  for (size_t row = 1; row < peps0.lambda_vert.rows() - 1; row++) {
    for (size_t col = 0; col < peps0.lambda_vert.cols(); col++) {
      Tensor &lambda = peps0.lambda_vert({row, col});
      Index<QNT> up_index = InverseIndex(peps0.Gamma({row - 1, col}).GetIndex(1));
      Index<QNT> dn_index = InverseIndex(peps0.Gamma({row, col}).GetIndex(3));
      lambda = Tensor({up_index, dn_index});
      for (size_t i = 0; i < up_index.dim(); i++)
        lambda({i, i}) = 1.0;
    }
  }

  for (size_t row = 0; row < peps0.lambda_horiz.rows(); row++) {
    for (size_t col = 1; col < peps0.lambda_horiz.cols() - 1; col++) {
      Tensor &lambda = peps0.lambda_horiz({row, col});
      Index<QNT> left_index = peps0.Gamma({row, col}).GetIndex(0);
      Index<QNT> right_index = InverseIndex(left_index);
      lambda = Tensor({left_index, right_index});
      for (size_t i = 0; i < left_index.dim(); i++)
        lambda({i, i}) = 1.0;
    }
  }

  auto su_exe = new qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, QNT>(update_para, peps0,
                                                                               ham_hei_nn);
  su_exe->Execute();
  auto tps2 = qlpeps::TPS<TenElemT, QNT>(su_exe->GetPEPS());
  SplitIndexTPS<TenElemT, QNT> split_index_tps2(tps2);
  split_index_tps2.Dump();
  return 0;
}