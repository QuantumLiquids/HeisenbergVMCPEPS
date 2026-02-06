// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-22
*
* Description: Unified Simple Update entry (Square Heisenberg/XY, Triangle Heisenberg).
*/

//#define PLAIN_TRANSPOSE 1

#include "qlpeps/algorithm/simple_update/square_lattice_nn_simple_update.h"
#include "qlpeps/algorithm/simple_update/square_lattice_nnn_simple_update.h"
#include "qlpeps/algorithm/simple_update/triangle_nn_on_sqr_peps_simple_update.h"
#include "qlpeps/qlpeps.h"
#include "qlpeps/api/conversions.h"
#include "./qldouble.h"
#include "./common_params.h"
#include <fstream>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <physics_params.json> <simple_update_algorithm_params.json>" << std::endl;
    return -1;
  }
  
  heisenberg_params::SimpleUpdateParams params(argv[1], argv[2]);

  // Model/data dispatch hints (Square XY handled by separate gate below)
  const std::string model = params.physical_params.ModelType;
  const bool is_triangle = (model.find("Triangle") != std::string::npos);
  const bool is_xy = (!is_triangle && (model == "SquareXY"));
  const auto bc = params.physical_params.BoundaryCondition;
  if (is_triangle && bc == qlpeps::BoundaryCondition::Periodic) {
    std::cerr << "ERROR: Triangle simple update with PBC is not supported." << std::endl;
    return -2;
  }

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

  // Heisenberg NN two-site gate
  ham_hei_nn({0, 0, 0, 0}) = 0.25;
  ham_hei_nn({1, 1, 1, 1}) = 0.25;
  ham_hei_nn({1, 1, 0, 0}) = -0.25;
  ham_hei_nn({0, 0, 1, 1}) = -0.25;
  ham_hei_nn({0, 1, 1, 0}) = 0.5;
  ham_hei_nn({1, 0, 0, 1}) = 0.5;

  // XY gate (SxSx + SySy) only has flip terms; no SzSz diagonal contribution
  Tensor ham_xy_nn = Tensor({pb_in, pb_out, pb_in, pb_out});
  ham_xy_nn({0, 1, 1, 0}) = 0.5; // |01><10|
  ham_xy_nn({1, 0, 0, 1}) = 0.5; // |10><01|

  // Select NN two-site Hamiltonian by ModelType
  Tensor ham_nn = is_xy ? ham_xy_nn : ham_hei_nn;
  auto ham_nnn = params.physical_params.J2 * ham_nn;

  // Triangle requires additional three-site term
  Tensor ham_hei_tri;
  if (is_triangle) {
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
    ham_hei_tri = ham_hei_tri_terms[0] + ham_hei_tri_terms[1] + ham_hei_tri_terms[2];
  }

  qlten::hp_numeric::SetTensorManipulationThreads(params.numerical_params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.numerical_params.Dmin, params.numerical_params.Dmax,
                                       params.numerical_params.TruncErr);

  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, params.physical_params.Ly, params.physical_params.Lx, bc);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    if (is_triangle) {
      // Three-sublattice ordered initial state for triangle lattice
      for (size_t y = 0; y < params.physical_params.Ly; y++) {
        for (size_t x = 0; x < params.physical_params.Lx; x++) {
          size_t sublattice_num;
          if (y >= x) {
            sublattice_num = (y - x) % 3;
          } else {
            sublattice_num = (3 - (x - y) % 3) % 3;
          }
          switch (sublattice_num) {
            case 0: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = 0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = 1; // spin down
                    break;
            case 1: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = -std::sqrt(3.0) / 2.0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
                    break;
            case 2: peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = std::sqrt(3.0) / 2.0;
                    peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
                    break;
          }
        }
      }
    } else {
      std::vector<std::vector<size_t>> activates(params.physical_params.Ly, std::vector<size_t>(params.physical_params.Lx));
      for (size_t y = 0; y < params.physical_params.Ly; y++) {
        for (size_t x = 0; x < params.physical_params.Lx; x++) {
          size_t sz_int = x + y;
          activates[y][x] = sz_int % 2;
        }
      }
      peps0.Initial(activates);
    }
  }

  std::unique_ptr<qlpeps::SimpleUpdateExecutor<TenElemT, QNT>> su_exe;
  if (is_triangle) {
    su_exe = std::make_unique<qlpeps::TriangleNNModelSquarePEPSSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_hei_nn, ham_hei_tri);
  } else if (std::abs(params.physical_params.J2) < 1e-15) {
    su_exe = std::make_unique<qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_nn);
  } else {
    su_exe = std::make_unique<qlpeps::SquareLatticeNNNSimpleUpdateExecutor<TenElemT, QNT>>(update_para, peps0, ham_nn, ham_nnn);
  }

  su_exe->Execute();
  
  // Convert to TPS (in-memory) and normalize, then generate SplitIndexTPS and save (TPS file not stored)
  auto peps = su_exe->GetPEPS();
  auto tps = qlpeps::ToTPS<TenElemT, QNT>(peps);
  for (auto &tensor : tps) {
      tensor.Normalize();
  }
  
  // Convert TPS to SplitIndexTPS and save for VMC use under final/lowest convention
  qlpeps::SplitIndexTPS<TenElemT, QNT> sitps = qlpeps::ToSplitIndexTPS<TenElemT, QNT>(tps);
  std::string sitps_final = "tpsfinal";
  sitps.Dump(sitps_final);

  // Save PEPS format for compatibility
  su_exe->DumpResult(peps_path, true);
  
  std::cout << "Simple Update completed." << std::endl;
  std::cout << "SplitIndexTPS saved to: " << sitps_final << std::endl;
  std::cout << "PEPS saved to: " << peps_path << std::endl;

  return 0;
}
