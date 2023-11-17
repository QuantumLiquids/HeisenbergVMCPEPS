//
// Created by haoxinwang on 02/11/2023.
//

#ifndef HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MEASUREMENT_SOLVER_H
#define HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MEASUREMENT_SOLVER_H


#include "gqten/gqten.h"
#include "gqpeps/two_dim_tn/tps/split_index_tps.h"      //SplitIndexTPS
#include "gqpeps/algorithm/vmc_update/tps_sample.h"     //TPSSample

namespace gqpeps {


//template<typename TenElemT, typename QNT>
//class MeasurementSolver {
//  using SITPS = SplitIndexTPS<TenElemT, QNT>;
// public:
//  MeasurementSolver(void) = default;
//
//  virtual TenElemT SampleMeasure(
//      const SITPS *sitps,
//      TPSSample<TenElemT, QNT> *tps_sample,
//      TensorNetwork2D<TenElemT, QNT> &hole_res  // the return value
//  ) {
//    TenElemT energy(0);
//    return energy;
//  }
//
// protected:
//};

template<typename TenElemT, typename QNT>
class KagomeSpinOneHalfHeisenbergMeasurementSolver {
  using SITPS = SplitIndexTPS<TenElemT, QNT>;
 public:
  KagomeSpinOneHalfHeisenbergMeasurementSolver(void) = default;

  TenElemT SampleMeasure(
      const SITPS *sitps,
      TPSSample<TenElemT, QNT> *tps_sample,
      std::vector<bool> &local_sz, //return 1
      std::vector<double> &energy_bond //return 2
  );
};


template<typename TenElemT, typename QNT>
TenElemT KagomeSpinOneHalfHeisenbergMeasurementSolver<TenElemT, QNT>::SampleMeasure(const SITPS *split_index_tps,
                                                                                    TPSSample<TenElemT, QNT> *tps_sample,
                                                                                    std::vector<bool> &local_sz,
                                                                                    std::vector<double> &energy_bond) {
  TenElemT energy(0);
  size_t ly = split_index_tps->rows(), lx = split_index_tps->cols();
  local_sz = std::vector<bool>(ly * lx * 3);
  energy_bond.clear();
  energy_bond.reserve(ly * lx * 6);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const TruncatePara &trunc_para = TPSSample<TenElemT, QNT>::trun_para;
  TenElemT inv_psi = 1.0 / (tps_sample->amplitude);
  tn.GenerateBMPSApproach(UP, trunc_para);
  for (size_t row = 0; row < tn.rows(); row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    for (size_t col = 0; col < tn.cols(); col++) {
      const SiteIdx site1 = {row, col};
      size_t config1 = config(site1);
      size_t site_1d_idx = (row * lx + col) * 3;
      local_sz[site_1d_idx] = config1 & 1; //left upper site
      local_sz[site_1d_idx + 1] = config1 >> 1 & 1;//lower site
      local_sz[site_1d_idx + 2] = config1 >> 2 & 1;//right site
      //      size_t config_left_upper = config(site1) % 2;
//      size_t config_lower = (config(site1) / 2) % 2;
//      size_t config_right = (config(site1) / 4);
      //Calculate the energy inner the three site of pseudo-site site1;
      if (config(site1) == 0 || config(site1) == 7) {
        energy += 0.75;
        energy_bond.push_back(0.25);
        energy_bond.push_back(0.25);
        energy_bond.push_back(0.25);
      } else {
        // 1->2->4->1
        // 3->6->5->3
        size_t rotate_config1 = config1 / 4 + 2 * (config1 % 4); //counter clockwise rotation
        size_t rotate_config2 = rotate_config1 / 4 + 2 * (rotate_config1 % 4);
        TenElemT psi_rotate1 = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[rotate_config1], HORIZONTAL);
        TenElemT psi_rotate2 = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[rotate_config2], HORIZONTAL);
        energy += -0.25 + (psi_rotate1 + psi_rotate2) * inv_psi * 0.5;
        if (config1 == 1 || config1 == 6) {
          energy_bond.push_back(-0.25 + psi_rotate1 * inv_psi * 0.5); //vertical bond in the triangular
          energy_bond.push_back(+0.25); //skew bond in the triangular
          energy_bond.push_back(-0.25 + psi_rotate2 * inv_psi * 0.5); //horizontal bond in the triangular
        } else if (config1 == 2 || config1 == 5) {
          energy_bond.push_back(-0.25 + psi_rotate2 * inv_psi * 0.5);
          energy_bond.push_back(-0.25 + psi_rotate1 * inv_psi * 0.5);
          energy_bond.push_back(+0.25);
        } else if (config1 == 4 || config1 == 3) {
          energy_bond.push_back(+0.25);
          energy_bond.push_back(-0.25 + psi_rotate2 * inv_psi * 0.5);
          energy_bond.push_back(-0.25 + psi_rotate1 * inv_psi * 0.5);
        }
      }

      if (col < tn.cols() - 1) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};
        size_t config2 = config(site2);
        if ((config1 >> 2 & 1) == (config2 & 1)) {
          energy += 0.25;
          energy_bond.push_back(0.25);
        } else {
          size_t ex_config1 = config1 ^ (1 << 2);
          size_t ex_config2 = config2 ^ 1;
          TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                  (*split_index_tps)(site1)[ex_config1],
                                                  (*split_index_tps)(site2)[ex_config2]);
          energy += (-0.25 + psi_ex * inv_psi * 0.5);
          energy_bond.push_back(-0.25 + psi_ex * inv_psi * 0.5);
        }
        tn.BTenMoveStep(RIGHT);
      }
    }
    if (row < tn.rows() - 1) {
      //Calculate diagonal energy contribution
      tn.InitBTen2(LEFT, row);
      tn.GrowFullBTen2(RIGHT, row, 2, true);

      for (size_t col = 0; col < tn.cols() - 1; col++) {
        //Calculate diagonal energy contribution
        SiteIdx site1 = {row + 1, col}; //left-down
        SiteIdx site2 = {row, col + 1}; //right-up
        size_t config1 = config(site1);
        size_t config2 = config(site2);
        if ((config1 >> 2 & 1) == (config2 >> 1 & 1)) {
          energy += 0.25;
          energy_bond.push_back(0.25);
        } else {
          size_t ex_config1 = config1 ^ (1 << 2);
          size_t ex_config2 = config2 ^ (1 << 1);
          TenElemT psi_ex = tn.ReplaceNNNSiteTrace({row, col},
                                                   LEFTDOWN_TO_RIGHTUP,
                                                   HORIZONTAL,
                                                   (*split_index_tps)(site1)[ex_config1],  //the tensor at left
                                                   (*split_index_tps)(site2)[ex_config2]);
          energy += (-0.25 + psi_ex * inv_psi * 0.5);
          energy_bond.push_back(-0.25 + psi_ex * inv_psi * 0.5);
        }
        tn.BTen2MoveStep(RIGHT, row);
      }
      tn.BMPSMoveStep(DOWN, trunc_para);
    }
  }

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    for (size_t row = 0; row < tn.rows() - 1; row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};
      size_t config1 = config(site1);
      size_t config2 = config(site2);
      if ((config1 >> 1 & 1) == (config2 & 1)) {
        energy += 0.25;
        energy_bond.push_back(0.25);
      } else {
        size_t ex_config1 = config1 ^ 1 << 1;
        size_t ex_config2 = config2 ^ 1;
        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                (*split_index_tps)(site1)[ex_config1],
                                                (*split_index_tps)(site2)[ex_config2]);
        energy += (-0.25 + psi_ex * inv_psi * 0.5);
        energy_bond.push_back(-0.25 + psi_ex * inv_psi * 0.5);
      }
      if (row < tn.rows() - 2) {
        tn.BTenMoveStep(DOWN);
      }
    }
    if (col < tn.cols() - 1) {
      tn.BMPSMoveStep(RIGHT, trunc_para);
    }
  }
  if (energy < -1.0e8) {
    std::cout << "Warning: sample's energy = " << energy << std::endl;
  }
  return energy;
}


}//gqpeps
#endif //HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MEASUREMENT_SOLVER_H
