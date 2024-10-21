/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-28
*
* Description: Simple Update for Heisenberg model.
*/


#define PLAIN_TRANSPOSE 1

#include "qlpeps/algorithm/loop_update/loop_update.h"
#include "./qldouble.h"
#include "./params_parser.h"

using LoopGateT = qlpeps::LoopGates<Tensor>;

#ifdef U1SYM
IndexT vb_out = IndexT({QNSctT(U1QN(0), 2),
                        QNSctT(U1QN(-1), 1),
                        QNSctT(U1QN(1), 1)},
                       TenIndexDirType::OUT
);
#else
IndexT vb_out = IndexT({
                           QNSctT(U1QN(0), 4)},
                       qlten::TenIndexDirType::OUT
);
#endif
IndexT vb_in = InverseIndex(vb_out);

LoopGateT GenerateTriangleHeisenbergLoopGates(
    const double tau, // imaginary time
    const size_t n0, const size_t n1, const size_t n2, const size_t n3  //bond share number
    //n0 bond is the upper horizontal bond in the loop
) {
  const std::vector<size_t> ns = {n0, n1, n2, n3};
  LoopGateT gates;

  for (size_t i = 0; i < 4; i++) {
    gates[i] = Tensor({vb_in, pb_in, pb_out, vb_out});
    Tensor &gate = gates[i];
    //Id
    gate({0, 0, 0, 0}) = 1.0;
    gate({0, 1, 1, 0}) = 1.0;
    //-s_z * tau
    gate({0, 0, 0, 1}) = -0.5 * tau / double(ns[i]);
    gate({0, 1, 1, 1}) = 0.5 * tau / double(ns[i]);
    //s_z
    gate({1, 0, 0, 0}) = 0.5;
    gate({1, 1, 1, 0}) = -0.5;

    //-s^+ * tau/2
    gate({0, 0, 1, 2}) = -1.0 * tau / double(ns[i]) / 2.0;
    //s^-
    gate({2, 1, 0, 0}) = 1.0;

    //-s^- * tau/2
    gate({0, 1, 0, 3}) = -1.0 * tau / double(ns[i]) / 2.0;
    //s^+
    gate({3, 0, 1, 0}) = 1.0;
  }

  for (auto i : {0}) {
    Tensor &gate = gates[i];
    gate({1, 0, 0, 1}) = double(ns[3]);
    gate({1, 1, 1, 1}) = double(ns[3]);

    gate({2, 0, 0, 2}) = double(ns[3]);
    gate({2, 1, 1, 2}) = double(ns[3]);

    gate({3, 0, 0, 3}) = double(ns[3]);
    gate({3, 1, 1, 3}) = double(ns[3]);
  }
  return gates;
}

void GenerateTriangleHeisenbergAllEvolveGates(
    const double tau, // imaginary time
    qlpeps::DuoMatrix<LoopGateT> &evolve_gates //output
) {
  //corner
  size_t Ly = evolve_gates.rows() + 1;
  size_t Lx = evolve_gates.cols() + 1; //physical size
  evolve_gates({0, 0}) = GenerateTriangleHeisenbergLoopGates(tau,
                                                             1, 2, 2, 1);
  evolve_gates({0, Lx - 2}) = GenerateTriangleHeisenbergLoopGates(tau,
                                                                  1, 1, 2, 2);
  evolve_gates({Ly - 2, 0}) = GenerateTriangleHeisenbergLoopGates(tau,
                                                                  2, 2, 1, 1);
  evolve_gates({Ly - 2, Lx - 2}) = GenerateTriangleHeisenbergLoopGates(tau,
                                                                       2, 1, 1, 2);

  auto gates_upper = GenerateTriangleHeisenbergLoopGates(tau,
                                                         1, 2, 2, 2);
  auto gates_lower = GenerateTriangleHeisenbergLoopGates(tau,
                                                         2, 2, 1, 2);
  for (size_t col = 1; col < Lx - 2; col++) {
    evolve_gates({0, col}) = gates_upper;
    evolve_gates({Ly - 2, col}) = gates_lower;
  }

  auto gates_left = GenerateTriangleHeisenbergLoopGates(tau,
                                                        2, 2, 2, 1);
  auto gates_middle = GenerateTriangleHeisenbergLoopGates(tau,
                                                          2, 2, 2, 2);
  auto gates_right = GenerateTriangleHeisenbergLoopGates(tau,
                                                         2, 1, 2, 2);
  for (size_t row = 1; row < Ly - 2; row++) {
    evolve_gates({row, 0}) = gates_left;
    evolve_gates({row, Lx - 2}) = gates_right;
  }
  for (size_t col = 1; col < Lx - 2; col++) {
    for (size_t row = 1; row < Ly - 2; row++) {
      evolve_gates({row, col}) = gates_middle;
    }
  }
}

int main(int argc, char **argv) {
  using namespace qlpeps;
  SimpleUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(1);
  omp_set_num_threads(params.ThreadNum);
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;

  qlpeps::DuoMatrix<LoopGateT> evolve_gates(Ly - 1, Lx - 1);
  GenerateTriangleHeisenbergAllEvolveGates(params.Tau, evolve_gates);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);

  ArnoldiParams arnoldi_params(1e-10, 200);
  double fet_tol = 1e-12;
  double fet_max_iter = 30;
  ConjugateGradientParams cg_params(100, 1e-10, 20, 0.0);

  FullEnvironmentTruncateParams fet_params(params.Dmin, params.Dmax, params.TruncErr,
                                           fet_tol, fet_max_iter,
                                           cg_params);

  qlpeps::SquareLatticePEPS<TenElemT, U1QN> peps0(pb_out, params.Ly, params.Lx);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    /**** Initialized as Three-Sublattice Order ****/
    for (size_t y = 0; y < params.Ly; y++) {
      for (size_t x = 0; x < params.Lx; x++) {
        size_t sublattice_num;
        if (y >= x) {
          sublattice_num = (y - x) % 3;
        } else {
          sublattice_num = (3 - (x - y) % 3) % 3;
        }
        switch (sublattice_num) {
          case 0:peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = 0;
            peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = 1; // spin down;
            break;
          case 1:peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = -std::sqrt(3.0) / 2.0;
            peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
            break;
          case 2:peps0.Gamma({y, x})({0, 0, 0, 0, 0}) = std::sqrt(3.0) / 2.0;
            peps0.Gamma({y, x})({0, 0, 0, 0, 1}) = -1.0 / 2.0;
            break;
        }
      }
    }
  }
  auto *loop_exe = new LoopUpdateExecutor<QLTEN_Double, U1QN>(LoopUpdateTruncatePara(
                                                                  arnoldi_params,
                                                                  1e-8,
                                                                  fet_params),
                                                              params.Step,
                                                              params.Tau,
                                                              evolve_gates,
                                                              peps0);

  loop_exe->Execute();
  auto peps_res = loop_exe->GetPEPS();
  loop_exe->DumpResult(peps_path, true);
  auto tps = qlpeps::TPS<TenElemT, U1QN>(peps_res);
  tps.Dump();
  delete loop_exe;
  return 0;
}