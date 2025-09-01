// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang <wanghaoxin1996@gmail.com>
* Creation Date: 2025-09-01
*
* Description: Convert legacy TPS data to SplitIndexTPS and dump to directory.
* Usage:
*   ./convert_tps_to_sitps <physics_params.json> [output_dir]
* Notes:
*   - Requires building without U1 symmetry for square-lattice Heisenberg example.
*   - When output_dir is omitted, defaults to "tpsfinal".
*/

#include "./qldouble.h"
#include "./params_parser.h"
#include "qlpeps/qlpeps.h"
#include "qlpeps/two_dim_tn/tps/tps.h"
#include "qlpeps/two_dim_tn/tps/split_index_tps.h"
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cout << "Usage: " << argv[0] << " <physics_params.json> [output_dir]" << std::endl;
    return -1;
  }

  const char *physics_file = argv[1];
  std::string output_dir = (argc == 3) ? std::string(argv[2]) : std::string("tpsfinal");
  if (!output_dir.empty() && output_dir.back() == '/') output_dir.pop_back();

  heisenberg_params::PhysicalParams phys(physics_file);

  try {
    qlten::hp_numeric::SetTensorManipulationThreads(1);

    // Load legacy TPS from current working directory using standard naming
    qlpeps::TPS<TenElemT, QNT> tps(phys.Ly, phys.Lx);
    bool ok = tps.Load();
    if (!ok) {
      std::cerr << "ERROR: Failed to load TPS files in current directory."
                << " Expect files like tps_*.qlten present." << std::endl;
      return -2;
    }

    // Normalize site tensors (safety)
    for (auto &tensor : tps) {
      tensor.Normalize();
    }

    // Convert to SplitIndexTPS and dump
    qlpeps::SplitIndexTPS<TenElemT, QNT> sitps(tps);
    sitps.Dump(output_dir);

    std::cout << "Converted TPS (" << phys.Ly << "x" << phys.Lx
              << ") to SplitIndexTPS at '" << output_dir << "'." << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return -3;
  }

  return 0;
}


