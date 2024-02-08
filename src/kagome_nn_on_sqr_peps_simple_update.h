// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-28
*
* Description: GraceQ/VMC-SquareLatticePEPS project.
*              Simple Update for nearest-neighbor interaction kagome lattice models in square lattice PEPS
*/



#ifndef GRACEQ_VMC_PEPS_Kagome_NN_ON_SQR_PEPS_SIMPLE_UPDATE_H
#define GRACEQ_VMC_PEPS_Kagome_NN_ON_SQR_PEPS_SIMPLE_UPDATE_H

#include "qlpeps/algorithm/simple_update/simple_update.h"

namespace qlpeps {

using namespace qlten;


template<typename TenElemT, typename QNT>
class KagomeNNModelSquarePEPSSimpleUpdateExecutor : public SimpleUpdateExecutor<TenElemT, QNT> {
  using Tensor = QLTensor<TenElemT, QNT>;
  using PEPST = SquareLatticePEPS<TenElemT, QNT>;
 public:
  /**
   * @param ham_nn     nearest-neighbor interaction, only consider top and left edge on the lattice
   * @param ham_tri    three-site triangle interaction term
   */
  KagomeNNModelSquarePEPSSimpleUpdateExecutor(const SimpleUpdatePara &update_para,
                                              const PEPST &peps_initial,
                                              const Tensor &ham_nn,
                                              const Tensor &ham_tri,
                                              const bool remove_corner) :
      SimpleUpdateExecutor<TenElemT, QNT>(update_para, peps_initial), \
         remove_corner_(remove_corner), ham_nn_(ham_nn), ham_tri_(ham_tri) {}

 private:
  void SetEvolveGate_(void) override {
    evolve_gate_nn_ = TaylorExpMatrix(this->update_para.tau, ham_nn_);
    evolve_gate_tri_ = TaylorExpMatrix(this->update_para.tau, ham_tri_);
  }

  double SimpleUpdateSweep_(void) override;

  bool remove_corner_;  // geometry feature
  Tensor ham_nn_;
  Tensor ham_tri_;
  Tensor evolve_gate_nn_;
  Tensor evolve_gate_tri_;
};

///< note PEPS linear size = 2 * linear size of system counting in unit cell
template<typename TenElemT, typename QNT>
double KagomeNNModelSquarePEPSSimpleUpdateExecutor<TenElemT, QNT>::SimpleUpdateSweep_(void) {
  Timer simple_update_sweep_timer("simple_update_sweep");
  SimpleUpdateTruncatePara para(this->update_para.Dmin, this->update_para.Dmax, this->update_para.Trunc_err);
  double norm = 1.0;
  double e0 = 0.0;

  //first row, horizontal bond link
  for (size_t col = 1; col < this->lx_ - 2; col += 2) {
    norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_, {0, col}, HORIZONTAL, para);
    e0 += -std::log(norm) / this->update_para.tau;
  }

  size_t upper_left_triangle_row_bound, upper_left_triangle_col_bound;
  if (remove_corner_) {
    upper_left_triangle_row_bound = this->ly_ - 3;
    upper_left_triangle_col_bound = this->lx_ - 3;
  } else {
    upper_left_triangle_row_bound = this->ly_ - 1;
    upper_left_triangle_col_bound = this->lx_ - 1;
  }
  for (size_t row = 0; row < upper_left_triangle_row_bound; row += 2) {
    for (size_t col = 0; col < upper_left_triangle_col_bound; col += 2) {
      norm = this->peps_.UpperLeftTriangleProject(evolve_gate_tri_, {row, col}, para);
      e0 += -std::log(norm) / this->update_para.tau;
    }
  }

  //first column, vertical bond link
  for (size_t row = 1; row < this->ly_ - 2; row += 2) {
    norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_, {row, 0}, VERTICAL, para);
    e0 += -std::log(norm) / this->update_para.tau;
  }

  for (size_t row = 1; row < this->ly_ - 2; row += 2) {
    for (size_t col = 2; col < this->lx_ - 1; col += 2) {
      norm = this->peps_.LowerRightTriangleProject(evolve_gate_tri_, {row, col}, para);
      e0 += -std::log(norm) / this->update_para.tau;
    }
  }
  if (remove_corner_) {
    //last row, horizontal bond link
    for (size_t col = 0; col < this->lx_ - 3; col += 2) {
      norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_, {this->ly_ - 2, col}, HORIZONTAL, para);
      e0 += -std::log(norm) / this->update_para.tau;
    }
    //last column,vertical bond link
    for (size_t row = 0; row < this->ly_ - 3; row += 2) {
      norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_, {row, this->lx_ - 2}, VERTICAL, para);
      e0 += -std::log(norm) / this->update_para.tau;
    }
  }

  double sweep_time = simple_update_sweep_timer.Elapsed();
  auto [dmin, dmax] = this->peps_.GetMinMaxBondDim();
  std::cout << "Estimated E0 =" << std::setw(15) << std::setprecision(kEnergyOutputPrecision) << std::fixed
            << std::right << e0
            << " Dmin/Dmax = " << std::setw(2) << std::right << dmin << "/" << std::setw(2) << std::left << dmax
            << " SweepTime = " << std::setw(8) << sweep_time
            << std::endl;
  return norm;
}
}


#endif//GRACEQ_VMC_PEPS_Kagome_NN_ON_SQR_PEPS_SIMPLE_UPDATE_H

