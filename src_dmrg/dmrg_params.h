// SPDX-License-Identifier: LGPL-3.0-only

/**
 * @file dmrg_params.h
 * @brief DMRG-specific case parameter parser.
 */

#ifndef HEISENBERGVMCPEPS_DMRG_PARAMS_H
#define HEISENBERGVMCPEPS_DMRG_PARAMS_H

#include "qlmps/qlmps.h"

#include <cstddef>
#include <string>
#include <vector>

/**
 * @brief Legacy-compatible DMRG parameter parser.
 *
 * Required keys in `CaseParams`:
 * - `Ly`, `Lx`, `RemoveCorner`, `J2`, `Sweeps`, `Dmin`, `Dmax`
 * - `CutOff`, `LanczErr`, `MaxLanczIter`, `Threads`, `noise`
 */
struct DMRGCaseParams : public qlmps::CaseParamsParserBasic {
  explicit DMRGCaseParams(const char *params_file)
      : qlmps::CaseParamsParserBasic(params_file),
        Geometry(),
        Ly(0),
        Lx(0),
        RemoveCorner(false),
        J1(0.0),
        J2(0.0),
        J3(0.0),
        Dzz(0.0),
        Sweeps(0),
        Dmin(0),
        Dmax(0),
        CutOff(0.0),
        LanczErr(0.0),
        MaxLanczIter(0),
        tau(0.0),
        M(0),
        ConvergeTolerance(0.0),
        TaylorOrder(0),
        TaylorErr(0.0),
        Threads(1),
        Perturbation(0.0),
        wavelength(0),
        noise(),
        SymmetryMode(0) {
    Ly = static_cast<size_t>(ParseInt("Ly"));
    Lx = static_cast<size_t>(ParseInt("Lx"));
    RemoveCorner = ParseBool("RemoveCorner");
    J2 = ParseDouble("J2");
    Sweeps = static_cast<size_t>(ParseInt("Sweeps"));
    Dmin = static_cast<size_t>(ParseInt("Dmin"));
    Dmax = static_cast<size_t>(ParseInt("Dmax"));
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = static_cast<size_t>(ParseInt("MaxLanczIter"));
    Threads = static_cast<size_t>(ParseInt("Threads"));
    noise = ParseDoubleVec("noise");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  bool RemoveCorner;
  double J1;
  double J2;
  double J3;
  double Dzz;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  double tau;
  size_t M;
  double ConvergeTolerance;
  size_t TaylorOrder;
  double TaylorErr;
  size_t Threads;
  double Perturbation;
  size_t wavelength;
  std::vector<double> noise;
  size_t SymmetryMode;  ///< Reserved legacy field.
};

#endif  // HEISENBERGVMCPEPS_DMRG_PARAMS_H
