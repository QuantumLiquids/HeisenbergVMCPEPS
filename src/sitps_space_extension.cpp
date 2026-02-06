//
// Created by haoxinwang on 02/01/2024.
// Duplicate the tensors in finite-size PEPS to extend the size of PEPS along one direction.
// Usage: ./extend <vmc params json file> <r or c> <start line> <increase width>";
// r or c means row or column, <vmc params json file> should give original size of PEPS,
//
// This program especially work for tensor without U1 case because for U1 it may induce the problem of indices dismatch.

#include "./qldouble.h"
#include "qlpeps/two_dim_tn/tps/split_index_tps.h"
#include "common_params.h"

using namespace qlpeps;

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cout << "Usage: ./extend <physics params json file> <r or c> <start line> <increase width>\n";
    return 1;
  }

  heisenberg_params::PhysicalParams physical_params(argv[1]);

  char direction = *argv[2];
  size_t start_line = std::atoi(argv[3]);
  size_t increase_width = std::atoi(argv[4]);

  SplitIndexTPS<TenElemT, QNT> original_tps(physical_params.Ly, physical_params.Lx);
  if (original_tps.Load()) {
    std::cout << "Loaded Original TPS." << std::endl;
  } else {
    std::cout << "Loading Original TPS fails!" << std::endl;
    exit(1);
  }

  if (direction == 'c') {
    SplitIndexTPS<TenElemT, QNT> extended_tps(physical_params.Ly, physical_params.Lx + increase_width);
    for (size_t row = 0; row < physical_params.Ly; row++) {
      for (size_t col = 0; col < start_line; col++) {
        extended_tps({row, col}) = original_tps({row, col});
      }
    }

    for (size_t row = 0; row < physical_params.Ly; row++) {
      for (size_t col = start_line; col < physical_params.Lx + increase_width; col++) {
        extended_tps({row, col}) = original_tps({row, col - increase_width});
      }
    }
    extended_tps.Dump();
  } else if (direction == 'r') {
    SplitIndexTPS<TenElemT, QNT> extended_tps(physical_params.Ly + increase_width, physical_params.Lx);
    for (size_t row = 0; row < start_line; row++) {
      for (size_t col = 0; col < physical_params.Lx; col++) {
        extended_tps({row, col}) = original_tps({row, col});
      }
    }

    for (size_t row = start_line; row < physical_params.Ly + increase_width; row++) {
      for (size_t col = 0; col < physical_params.Lx; col++) {
        extended_tps({row, col}) = original_tps({row - increase_width, col});
      }
    }
    extended_tps.Dump();
  } else {
    std::cout << "Invalid direction parameter. Should be 'r' or 'c'.\n";
    return 1;
  }

  return 0;
}