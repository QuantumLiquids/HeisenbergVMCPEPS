//
// Created by haoxinwang on 02/01/2024.
//


#include "./qldouble.h"
#include "qlpeps/two_dim_tn/tps/split_index_tps.h"
#include "params_parser.h"
#include "myutil.h"

using namespace qlpeps;

int main(int argc, char **argv) {

  VMCUpdateParams params(argv[1]);

  SplitIndexTPS < TenElemT, U1QN > original_tps(params.Ly, params.Lx);
  original_tps.Load();

  SplitIndexTPS < TenElemT, U1QN > extended_tps(params.Ly, params.Lx + 3);
  for (size_t row = 0; row < params.Ly; row++) {
    for (size_t col = 0; col < 4; col++) {
      extended_tps({row, col}) = original_tps({row, col});
    }
  }

  for (size_t row = 0; row < params.Ly; row++) {
    for (size_t col = 4; col < params.Lx + 3; col++) {
      extended_tps({row, col}) = original_tps({row, col - 3});
    }
  }
  extended_tps.Dump();

  return 0;
}