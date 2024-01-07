/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-12-25
*
* Description: GraceQ/VMC-PEPS project. Model Energy Solver for spin-1/2 AFM Heisenberg model in square lattice
*/

#ifndef GRACEQ_VMC_PEPS_SPIN_ONEHALF_HEISENBERG_SQUARE_BOOST_H
#define GRACEQ_VMC_PEPS_SPIN_ONEHALF_HEISENBERG_SQUARE_BOOST_H

#include "gqpeps/algorithm/vmc_update/model_energy_solver.h"    //ModelEnergySolver

namespace gqpeps {
using namespace gqten;

size_t NumOfSpinDownInEvenSublattice(const Configuration &config) {
  size_t n;
  for (size_t row = 0; row < config.rows(); row += 2) {
    for (size_t col = 0; col < config.cols(); col += 2) {
      n += (1 - config({row, col}));  //suppose 0 means spin down;
    }
  }
  return n;
}

DuoMatrix<double> WeightAttachedToSite(const DuoMatrix<double> &u, const Configuration &config) {
  // total weights = sum of weight.
  // total weights is the exponent in exp
  // store weight relevant to site so that it's easier to get the total weight when exchange NN spin.
  DuoMatrix<double> weights(u.rows(), u.cols());
  for (size_t row1 = 0; row1 < u.rows(); row1++) {
    for (size_t col1 = 0; col1 < u.cols(); col1++) {
      const SiteIdx site1 = {row1, col1};
      const size_t sigma1 = config(site1);

      double weight = 0.0;
      for (size_t row2 = 0; row2 < u.rows(); row2++) {
        for (size_t col2 = 0; col2 < u.cols(); col2++) {
          const SiteIdx site2 = {row2, col2};
          if (site1 == SiteIdx(site2)) {
            continue;
          }
          weight += u(site2) * config(site2);
        }
      }
      weight *= (-0.25 * sigma1);
      weights(site1) = weight;
    }
  }
  return weights;
}

template<typename TenElemT, typename QNT>
class SpinOneHalfHeisenbergSquareBoost : public ModelEnergySolver<TenElemT, QNT> {
  using SITPS = SplitIndexTPS<TenElemT, QNT>;
 public:
  SpinOneHalfHeisenbergSquareBoost(void) = delete;

  SpinOneHalfHeisenbergSquareBoost(const size_t lx, const size_t ly)
      : lx_(lx), ly_(ly), k_(ly, lx), gamma_(ly, lx), u_(ly, lx) {
    const size_t N = lx * ly;
    DuoMatrix<double> form_factors(ly_, lx_); //((1+gamma)/(1-gamma))^1/2 - 1  of k
    for (size_t row = 0; row < ly_; row++) {
      double ky = (double) (row) / double(ly_) * (2 * M_PI) - M_PI;
      for (size_t col = 0; col < lx_; col++) {
        double kx = (double) (col) / double(lx_) * (2 * M_PI) - M_PI;
        k_({row, col}) = std::make_pair(kx, ky);
        double gamma = 0.5 * (std::cos(kx) + std::cos(ky));
        gamma_({row, col}) = gamma;
        form_factors({row, col}) = std::sqrt((1 + gamma) / (1 - gamma)) - 1;
      }
    }

    for (size_t ry = 0; ry < ly_; ry++) {
      for (size_t rx = 0; rx < lx_; rx++) {
        //(rx, ry) are distance from site i to site j
        std::complex<double> u_ij = 0.0;
        for (size_t ky_idx = 0; ky_idx < ly_; ky_idx++) {
          for (size_t kx_idx = 0; kx_idx < lx_; kx_idx++) {
            auto k = k_({ky_idx, kx_idx});
            double form_factor = form_factors({ky_idx, kx_idx});
            double phase_factor = k.first * rx + k.second * ry;
            u_ij += form_factor * std::exp(-std::complex<double>(0.0, phase_factor));
          }
        }
        u_ij /= double(N);
        u_({ry, rx}) = u_ij.real(); // inversion symmetry of gamma(k) makes sure the imag part is 0.
      }
    }
  }

  TenElemT CalEnergyAndHoles(
      const SITPS *sitps,
      TPSSample<TenElemT, QNT> *tps_sample,
      TensorNetwork2D<TenElemT, QNT> &hole_res
  ) override;

  TenElemT CalEnergy(
      const SITPS *sitps,
      TPSSample<TenElemT, QNT> *tps_sample
  ) override;
 private:
  const size_t lx_; //cols
  const size_t ly_; //rows

  //Phys. Rev. B 40, 11437(R)
  DuoMatrix<std::pair<double, double>> k_; //momenta in F.B.Z., (kx, ky)
  DuoMatrix<double> gamma_; // gamma of k
  DuoMatrix<double> u_; // u of distance. fix site i = (0,0) due to translation invariance.
};

template<typename TenElemT, typename QNT>
TenElemT SpinOneHalfHeisenbergSquareBoost<TenElemT, QNT>::CalEnergyAndHoles(const SITPS *split_index_tps,
                                                                            TPSSample<TenElemT, QNT> *tps_sample,
                                                                            TensorNetwork2D<TenElemT, QNT> &hole_res) {
  TenElemT energy(0);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const DuoMatrix<double> weights = WeightAttachedToSite(u_, config);
  double exponent_sum = 0.0;
  for (size_t ry = 0; ry < ly_; ry++) {
    for (size_t rx = 0; rx < lx_; rx++) {
      exponent_sum += weights({ry, rx});
    }
  }
  int marshall_sign;
  size_t n = NumOfSpinDownInEvenSublattice(config);
  if (n % 2 == 1) {
    marshall_sign = -1;
  } else {
    marshall_sign = 1;
  }
  TenElemT psi0 = std::exp(exponent_sum) * marshall_sign;
  const BMPSTruncatePara &trunc_para = TPSSample<TenElemT, QNT>::trun_para;
  TenElemT inv_psi = 1.0 / (tps_sample->amplitude + psi0);
  tn.GenerateBMPSApproach(UP, trunc_para);
  for (size_t row = 0; row < tn.rows(); row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / (tps_sample->amplitude + psi0);
    for (size_t col = 0; col < tn.cols(); col++) {
      const SiteIdx site1 = {row, col};
      //Calculate the holes
      hole_res(site1) = Dag(tn.PunchHole(site1, HORIZONTAL));
      if (col < tn.cols() - 1) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};
        if (config(site1) == config(site2)) {
          energy += 0.25;
        } else {
          TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
          TenElemT psi_prime_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                        (*split_index_tps)(site1)[config(site2)],
                                                        (*split_index_tps)(site2)[config(site1)]);
          energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
        }
        tn.BTenMoveStep(RIGHT);
      }
      //PBC link... first we only write the VMC part
      const SiteIdx site2 = {row, 0};
      if (config(site1) == config(site2)) {
        energy += 0.25;
      } else {
        TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
        TenElemT psi_prime_ex = 0;
        energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
      }
    }
    if (row < tn.rows() - 1) {
      tn.BMPSMoveStep(DOWN, trunc_para);
    }
  }

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / (tps_sample->amplitude + psi0);
    for (size_t row = 0; row < tn.rows() - 1; row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};
      if (config(site1) == config(site2)) {
        energy += 0.25;
      } else {
        TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
        TenElemT psi_prime_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                      (*split_index_tps)(site1)[config(site2)],
                                                      (*split_index_tps)(site2)[config(site1)]);
        energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
      }
      if (row < tn.rows() - 2) {
        tn.BTenMoveStep(DOWN);
      }
    }
    //PBC link, only include the VMC part.
    size_t row = tn.rows() - 1;
    const SiteIdx site1 = {row, col};
    const SiteIdx site2 = {0, col};
    if (config(site1) == config(site2)) {
      energy += 0.25;
    } else {
      TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
      TenElemT psi_prime_ex = 0;
      energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
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

template<typename TenElemT, typename QNT>
TenElemT SpinOneHalfHeisenbergSquareBoost<TenElemT, QNT>::CalEnergy(const SITPS *split_index_tps,
                                                                    TPSSample<TenElemT, QNT> *tps_sample) {
  TenElemT energy(0);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const DuoMatrix<double> weights = WeightAttachedToSite(u_, config);
  double exponent_sum = 0.0;
  for (size_t ry = 0; ry < ly_; ry++) {
    for (size_t rx = 0; rx < lx_; rx++) {
      exponent_sum += weights({ry, rx});
    }
  }
  int marshall_sign;
  size_t n = NumOfSpinDownInEvenSublattice(config);
  if (n % 2 == 1) {
    marshall_sign = -1;
  } else {
    marshall_sign = 1;
  }
  TenElemT psi0 = std::exp(exponent_sum) * marshall_sign;
  const BMPSTruncatePara &trunc_para = TPSSample<TenElemT, QNT>::trun_para;
  TenElemT inv_psi = 1.0 / (tps_sample->amplitude + psi0);
  tn.GenerateBMPSApproach(UP, trunc_para);
  for (size_t row = 0; row < tn.rows(); row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / (tps_sample->amplitude + psi0);
    for (size_t col = 0; col < tn.cols(); col++) {
      const SiteIdx site1 = {row, col};
      if (col < tn.cols() - 1) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};
        if (config(site1) == config(site2)) {
          energy += 0.25;
        } else {
          TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
          TenElemT psi_prime_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                        (*split_index_tps)(site1)[config(site2)],
                                                        (*split_index_tps)(site2)[config(site1)]);
          energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
        }
        tn.BTenMoveStep(RIGHT);
      }
      //PBC link... first we only write the VMC part
      const SiteIdx site2 = {row, 0};
      if (config(site1) == config(site2)) {
        energy += 0.25;
      } else {
        TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
        TenElemT psi_prime_ex = 0;
        energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
      }
    }
    if (row < tn.rows() - 1) {
      tn.BMPSMoveStep(DOWN, trunc_para);
    }
  }

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / (tps_sample->amplitude + psi0);
    for (size_t row = 0; row < tn.rows() - 1; row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};
      if (config(site1) == config(site2)) {
        energy += 0.25;
      } else {
        TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
        TenElemT psi_prime_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                      (*split_index_tps)(site1)[config(site2)],
                                                      (*split_index_tps)(site2)[config(site1)]);
        energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
      }
      if (row < tn.rows() - 2) {
        tn.BTenMoveStep(DOWN);
      }
    }
    //PBC link, only include the VMC part.
    size_t row = tn.rows() - 1;
    const SiteIdx site1 = {row, col};
    const SiteIdx site2 = {0, col};
    if (config(site1) == config(site2)) {
      energy += 0.25;
    } else {
      TenElemT psi0_ex = std::exp(exponent_sum - 2 * weights(site1) - 2 * weights(site2)) * marshall_sign * (-1);
      TenElemT psi_prime_ex = 0;
      energy += (-0.25 + (psi0_ex + psi_prime_ex) * inv_psi * 0.5);
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

#endif //GRACEQ_VMC_PEPS_SPIN_ONEHALF_HEISENBERG_SQUARE_BOOST_H
