//
// Created by haoxinwang on 11/10/2023.
//

/**
 * The configuration of the combined pseudo-site is a number from 0 to 7.
 * After rewriting the configuration in a binary number,
 * the lowest bit corresponds to the left-upper site of the combined pseudo-site,
 * the next bit is lower site,
 * and the highest bit corresponds to the right site.
 */

#ifndef HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MODEL_SQRPEPS_SOLVER_H
#define HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MODEL_SQRPEPS_SOLVER_H

#include "gqpeps/two_dim_tn/peps/square_lattice_peps.h"
#include "gqpeps/algorithm/vmc_update/model_energy_solver.h"    //ModelEnergySolver
#include "gqpeps/algorithm/vmc_update/model_measurement_solver.h"    //

namespace gqpeps {
using namespace gqten;

template<typename TenElemT, typename QNT>
SplitIndexTPS<TenElemT, QNT> KagomeSquarePEPSToSplitIndexTPS(
    const SquareLatticePEPS<TenElemT, QNT> &peps
) {
  TPS<TenElemT, QNT> tps = TPS<TenElemT, QNT>(peps);
  size_t rows = peps.Rows() / 2, cols = peps.Cols() / 2;
  size_t phy_dim = tps({0, 0}).GetIndex(4).dim(); //2
  size_t combined_phy_dim = phy_dim * phy_dim * phy_dim;
  SplitIndexTPS<TenElemT, QNT> split_idx_tps(rows, cols);
  using Tensor = GQTensor<TenElemT, QNT>;
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < cols; col++) {
      Tensor left_upper = tps({2 * row, 2 * col});
      Tensor right_upper = tps({2 * row, 2 * col + 1});
      Tensor left_lower = tps({2 * row + 1, 2 * col});
      right_upper.RemoveTrivialIndexes({1, 3});
      left_lower.RemoveTrivialIndexes({0, 2});

      Tensor tmp, combined_tps_ten;
      Contract(&left_upper, {1}, &left_lower, {1}, &tmp);
      Contract(&tmp, {1}, &right_upper, {0}, &combined_tps_ten);
      combined_tps_ten.Transpose({0, 3, 5, 1, 2, 4, 6});
      /**
       * combined_tps_ten up to now
       *           3
       *           |
       * 0--combined_tps_ten--2
       *           |
       *           1
       *  and physical indexes:
       *    4: left_upper site
       *    5: lower site
       *    6: right site
       */
      combined_tps_ten.Normalize();
      std::vector<Index<QNT>> indexes = combined_tps_ten.GetIndexes();
      std::vector<Index<QNT>> phy_indexes_out(indexes.begin() + 4, indexes.end());
      std::vector<Index<QNT>> phy_indexes_in(3);
      for (size_t i = 0; i < 3; i++) {
        phy_indexes_in[i] = InverseIndex(phy_indexes_out[i]);
      }
      split_idx_tps({row, col}) = std::vector<Tensor>(combined_phy_dim);
      for (size_t dim = 0; dim < combined_phy_dim; dim++) {
        Tensor proj_ten(phy_indexes_in);
        proj_ten({dim % phy_dim, (dim / phy_dim) % (phy_dim), (dim / phy_dim / phy_dim) % (phy_dim)}) = 1.0;
        Contract(&combined_tps_ten, {4, 5, 6}, &proj_ten, {0, 1, 2}, &split_idx_tps({row, col})[dim]);
        if (split_idx_tps({row, col})[dim].GetQNBlkNum() == 0) {
          std::cout << "warning: Site (" << row << "," << col << "), dim = " << dim << " projected tensor is empty."
                    << std::endl;
        }
      }
    }
  }
  return split_idx_tps;
}

template<typename TenElemT, typename QNT>
class KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver : public ModelEnergySolver<TenElemT, QNT>,
                                                      ModelMeasurementSolver<TenElemT, QNT> {
  using SITPS = SplitIndexTPS<TenElemT, QNT>;
 public:
  using ModelEnergySolver<TenElemT, QNT>::ModelEnergySolver;

  KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver(bool remove_corner) : remove_corner_(remove_corner) {}

  template<typename WaveFunctionComponentType>
  ObservablesLocal<TenElemT> SampleMeasure(
      const SITPS *,
      WaveFunctionComponentType *
  );

  template<typename WaveFunctionComponentType, bool calchole = true>
  TenElemT CalEnergyAndHoles(
      const SITPS *sitps,
      WaveFunctionComponentType *tps_sample,
      TensorNetwork2D<TenElemT, QNT> &hole_res
  );

 private:
  bool remove_corner_ = true;
};

template<typename TenElemT, typename QNT>
template<typename WaveFunctionComponentType>
ObservablesLocal<TenElemT> KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver<TenElemT, QNT>::SampleMeasure(
    const SITPS *split_index_tps,
    WaveFunctionComponentType *tps_sample
) {
  ObservablesLocal<TenElemT> res;
  TenElemT &energy = res.energy_loc;
  std::vector<TenElemT> &local_sz = res.one_point_functions_loc;
  std::vector<TenElemT> &energy_bond = res.bond_energys_loc;
  size_t ly = split_index_tps->rows(), lx = split_index_tps->cols();
  local_sz = std::vector<TenElemT>(ly * lx * 3 - lx - ly);
  energy_bond.clear();
  energy_bond.reserve(ly * lx * 6);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const BMPSTruncatePara &trunc_para = tps_sample->trun_para;
  TenElemT inv_psi;
  tn.GenerateBMPSApproach(UP, trunc_para);
  for (size_t row = 0; row < tn.rows() - 1; row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    for (size_t col = 0; col < tn.cols() - 1; col++) {
      const SiteIdx site1 = {row, col};
      size_t config1 = config(site1);
      size_t site_1d_idx = row * (3 * lx - 1) + col * 3;
      local_sz[site_1d_idx] = config1 & 1; //left upper site
      local_sz[site_1d_idx + 1] = config1 >> 1 & 1;//lower site
      local_sz[site_1d_idx + 2] = config1 >> 2 & 1;//right site
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
    //last vertical bond in the incomplete triangle
    const SiteIdx site1 = {row, lx - 1};
    size_t config1 = config(site1);
    size_t site_1d_idx = (row * lx + lx - 1) * 3 - row;
    local_sz[site_1d_idx] = config1 & 1; //left upper site
    local_sz[site_1d_idx + 1] = config1 >> 1 & 1;//lower site

    bool site_upper = (config1 & 1);
    bool site_lower = (config1 >> 1 & 1);
    if (site_upper == site_lower) {
      energy += 0.25;
      energy_bond.push_back(0.25);
    } else {
      size_t ex_config = 2 * site_upper + site_lower + 4 * (config1 >> 2 & 1);
      TenElemT psi_ex = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[ex_config], HORIZONTAL);
      energy += (-0.25 + psi_ex * inv_psi * 0.5);
      energy_bond.push_back(-0.25 + psi_ex * inv_psi * 0.5);
    }

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
  //energy in last row
  size_t last_row_idx = ly - 1;
  tn.InitBTen(LEFT, last_row_idx);
  tn.GrowFullBTen(RIGHT, last_row_idx, 1, true);
  for (size_t col = 0; col < tn.cols() - 1; col++) {
    // bond energy inner PEPS site
    const SiteIdx site1 = {last_row_idx, col};
    size_t config1 = config(site1);
    size_t site_1d_idx = last_row_idx * (lx * 3 - 1) + col * 2;
    local_sz[site_1d_idx] = config1 & 1; //left upper site
    local_sz[site_1d_idx + 1] = config1 >> 2 & 1;//right site
    bool site_left = (config1 & 1);
    bool site_right = (config1 >> 2 & 1);
    if (site_left == site_right) {
      energy += 0.25;
      energy_bond.push_back(0.25);
    } else {
      size_t ex_config = 4 * site_left + 2 * (config1 >> 1 & 1) + site_right;
      TenElemT psi_ex = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[ex_config], HORIZONTAL);
      energy += (-0.25 + psi_ex * inv_psi * 0.5);
      energy_bond.push_back(-0.25 + psi_ex * inv_psi * 0.5);
    }

    //bond energy intra PEPS site
    const SiteIdx site2 = {last_row_idx, col + 1};
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
  size_t last_site_idx = last_row_idx * (lx * 3 - 1) + (lx - 1) * 2;
  local_sz[last_site_idx] = config({ly - 1, lx - 1}) & 1;

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / tps_sample->amplitude;
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
  for (auto &sz : local_sz) {
    sz -= 0.5;
  }
  if (energy < -1.0e8) {
    std::cout << "Warning: sample's energy = " << energy << std::endl;
  }
  return res;
}

template<typename TenElemT, typename QNT>
template<typename WaveFunctionComponentType, bool calchole>
TenElemT KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver<TenElemT, QNT>::CalEnergyAndHoles(const SITPS *split_index_tps,
                                                                                         WaveFunctionComponentType *tps_sample,
                                                                                         TensorNetwork2D<TenElemT,
                                                                                                         QNT> &hole_res) {
  TenElemT energy(0);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const BMPSTruncatePara &trunc_para = tps_sample->trun_para;
  TenElemT inv_psi;
  const size_t tri_right_bound = tn.cols() - (size_t) remove_corner_;
  const size_t tri_lower_bound = tn.rows() - (size_t) remove_corner_;
  tn.GenerateBMPSApproach(UP, trunc_para);
  for (size_t row = 0; row < tn.rows(); row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    for (size_t col = 0; col < tn.cols(); col++) {
      const SiteIdx site1 = {row, col};
      //Calculate the holes
      if constexpr (calchole)
        hole_res(site1) = Dag(tn.PunchHole(site1, HORIZONTAL));
      size_t config1 = config(site1);
      //      size_t config_left_upper = config(site1) % 2;
//      size_t config_lower = (config(site1) / 2) % 2;
//      size_t config_right = (config(site1) / 4);
      //Calculate the energy inner the three site of pseudo-site site1;
      if (col < tri_right_bound && row < tri_lower_bound) {
        if (config(site1) == 0 || config(site1) == 7) {
          energy += 0.75;
        } else {
          // 1->2->4->1
          // 3->6->5->3
          size_t rotate_config1 = config1 / 4 + 2 * (config1 % 4);
          size_t rotate_config2 = rotate_config1 / 4 + 2 * (rotate_config1 % 4);
          TenElemT psi_rotate1 = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[rotate_config1], HORIZONTAL);
          TenElemT psi_rotate2 = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[rotate_config2], HORIZONTAL);
          energy += -0.25 + (psi_rotate1 + psi_rotate2) * inv_psi * 0.5;
        }
      }

      if (remove_corner_ && col == tri_right_bound && row < tri_lower_bound) {
        bool site_upper = (config1 & 1);
        bool site_lower = (config1 >> 1 & 1);
        if (site_upper == site_lower) {
          energy += 0.25;
        } else {
          size_t ex_config = 2 * site_upper + site_lower + (config1 & 4);
          TenElemT psi_ex = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[ex_config], HORIZONTAL);
          energy += (-0.25 + psi_ex * inv_psi * 0.5);
        }
      }

      if (col < tn.cols() - 1) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};
        size_t config2 = config(site2);
        if ((config1 >> 2 & 1) == (config2 & 1)) {
          energy += 0.25;
        } else {
          size_t ex_config1 = config1 ^ (1 << 2);
          size_t ex_config2 = config2 ^ 1;
          TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                  (*split_index_tps)(site1)[ex_config1],
                                                  (*split_index_tps)(site2)[ex_config2]);
          energy += (-0.25 + psi_ex * inv_psi * 0.5);
        }

        if (remove_corner_ && row == tri_lower_bound) {
          bool site_left = (config1 & 1);
          bool site_right = (config1 >> 2 & 1);
          if (site_left == site_right) {
            energy += 0.25;
          } else {
            size_t ex_config = 4 * site_left + (config1 & 2) + site_right;
            TenElemT psi_ex = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[ex_config], HORIZONTAL);
            energy += (-0.25 + psi_ex * inv_psi * 0.5);
          }
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
        } else {
          size_t ex_config1 = config1 ^ (1 << 2);
          size_t ex_config2 = config2 ^ (1 << 1);
          TenElemT psi_ex = tn.ReplaceNNNSiteTrace({row, col},
                                                   LEFTDOWN_TO_RIGHTUP,
                                                   HORIZONTAL,
                                                   (*split_index_tps)(site1)[ex_config1],  //the tensor at left
                                                   (*split_index_tps)(site2)[ex_config2]);
          energy += (-0.25 + psi_ex * inv_psi * 0.5);
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
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    for (size_t row = 0; row < tn.rows() - 1; row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};
      size_t config1 = config(site1);
      size_t config2 = config(site2);
      if ((config1 >> 1 & 1) == (config2 & 1)) {
        energy += 0.25;
      } else {
        size_t ex_config1 = config1 ^ 1 << 1;
        size_t ex_config2 = config2 ^ 1;
        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                (*split_index_tps)(site1)[ex_config1],
                                                (*split_index_tps)(site2)[ex_config2]);
        energy += (-0.25 + psi_ex * inv_psi * 0.5);
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


#endif //HEISENBERGVMCPEPS_SPIN_ONEHALF_HEISENBERG_KAGOME_MODEL_SQRPEPS_SOLVER_H
