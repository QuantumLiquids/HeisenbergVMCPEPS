// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-12-19
*
* Description: Common physical and numerical parameters shared across different algorithms.
*/

#ifndef HEISENBERGVMCPEPS_COMMON_PARAMS_H
#define HEISENBERGVMCPEPS_COMMON_PARAMS_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/ond_dim_tn/boundary_mps/bmps.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"

namespace heisenberg_params {

/**
 * @brief Common physical parameters for lattice models
 */
struct PhysicalParams : public qlmps::CaseParamsParserBasic {
  PhysicalParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    J2 = ParseDouble("J2");
    RemoveCorner = ParseBool("RemoveCorner");
    // Optional, with sensible defaults documented in tutorials
    try { ModelType = ParseStr("ModelType"); } catch (...) { ModelType = "SquareHeisenberg"; }
    try { MCRestrictU1 = ParseBool("MCRestrictU1"); } catch (...) { MCRestrictU1 = true; }
  }

  size_t Lx;          ///< Lattice size in x direction
  size_t Ly;          ///< Lattice size in y direction
  double J2;          ///< J1-J2 model parameter (J1 = 1.0 by default)
  bool RemoveCorner;  ///< Whether to remove corner sites (for specific geometries)
  std::string ModelType; ///< Physics model identifier (e.g., SquareHeisenberg, SquareXY, TriangleHeisenberg)
  bool MCRestrictU1;     ///< Whether MC sampler restricts U1 (true by default)
};

/**
 * @brief Common numerical parameters for tensor network calculations
 */
struct NumericalParams : public qlmps::CaseParamsParserBasic {
  NumericalParams(const char *f) : CaseParamsParserBasic(f) {
    TruncErr = ParseDouble("TruncErr");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    ThreadNum = ParseInt("ThreadNum");
  }

  double TruncErr;    ///< Truncation error threshold
  size_t Dmin;        ///< Minimum bond dimension
  size_t Dmax;        ///< Maximum bond dimension
  size_t ThreadNum;   ///< Number of threads for tensor operations
};

/**
 * @brief Monte Carlo specific parameters
 */
struct MonteCarloNumericalParams : public qlmps::CaseParamsParserBasic {
  MonteCarloNumericalParams(const char *f) : CaseParamsParserBasic(f) {
    MC_samples = ParseInt("MC_samples");
    WarmUp = ParseInt("WarmUp");
    MCLocalUpdateSweepsBetweenSample = ParseInt("MCLocalUpdateSweepsBetweenSample");
  }

  size_t MC_samples;                        ///< Number of Monte Carlo samples
  size_t WarmUp;                           ///< Number of warm-up sweeps
  size_t MCLocalUpdateSweepsBetweenSample; ///< Sweeps between successive samples
};

/**
 * @brief Boundary MPS specific parameters
 */
struct BMPSParams : public qlmps::CaseParamsParserBasic {
  BMPSParams(const char *f) : CaseParamsParserBasic(f) {
    // Boundary MPS bond dimensions
    Db_max = ParseInt("Dbmps_max");
    try { Db_min = ParseInt("Dbmps_min"); } catch (...) { Db_min = Db_max; }
    MPSCompressScheme = static_cast<qlpeps::CompressMPSScheme>(ParseInt("MPSCompressScheme"));
    // Numerical controls specific to BMPS
    TruncErr = ParseDouble("TruncErr");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Db_min;                              ///< Minimum boundary MPS bond dimension
  size_t Db_max;                              ///< Maximum boundary MPS bond dimension
  qlpeps::CompressMPSScheme MPSCompressScheme; ///< MPS compression scheme
  double TruncErr;                             ///< Truncation error threshold for boundary MPS
  size_t ThreadNum;                            ///< Number of threads for tensor operations

  /**
   * @brief Create BMPSTruncatePara from these parameters
   */
  qlpeps::BMPSTruncatePara CreateTruncatePara(double trunc_err) const {
    return qlpeps::BMPSTruncatePara(Db_min, Db_max, trunc_err, MPSCompressScheme,
                                   std::make_optional<double>(trunc_err),
                                   std::make_optional<size_t>(10));
  }

  /**
   * @brief Create BMPSTruncatePara using internal TruncErr
   */
  qlpeps::BMPSTruncatePara CreateTruncatePara() const {
    return qlpeps::BMPSTruncatePara(Db_min, Db_max, TruncErr, MPSCompressScheme,
                                   std::make_optional<double>(TruncErr),
                                   std::make_optional<size_t>(10));
  }
};

} // namespace heisenberg_params

#endif //HEISENBERGVMCPEPS_COMMON_PARAMS_H
