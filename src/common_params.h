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
#include "qlpeps/one_dim_tn/boundary_mps/bmps.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "qlpeps/two_dim_tn/common/boundary_condition.h"
#include "qlpeps/two_dim_tn/tensor_network_2d/trg/trg_contractor.h"
#include <algorithm>
#include <cctype>
#include <optional>
#include <stdexcept>

namespace heisenberg_params {

inline std::string TrimAsciiWhitespace(std::string s) {
  auto is_space = [](unsigned char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
  };
  while (!s.empty() && is_space(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
  while (!s.empty() && is_space(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

inline qlpeps::BoundaryCondition ParseBoundaryCondition(const std::string &value) {
  std::string key = TrimAsciiWhitespace(value);
  std::transform(key.begin(), key.end(), key.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (key == "open" || key == "obc") return qlpeps::BoundaryCondition::Open;
  if (key == "periodic" || key == "pbc") return qlpeps::BoundaryCondition::Periodic;
  throw std::invalid_argument("BoundaryCondition must be Open/OBC or Periodic/PBC.");
}

inline qlpeps::CompressMPSScheme ParseCompressMPSScheme(const std::string &value) {
  std::string key = value;
  std::transform(key.begin(), key.end(), key.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (key == "svd" || key == "svd_compress" || key == "svdcompression" ||
      key == "svd compression") {
    return qlpeps::CompressMPSScheme::SVD_COMPRESS;
  }
  if (key == "var2" || key == "variational2site" || key == "variation2site" ||
      key == "two-site variational compression" || key == "two_site_variational") {
    return qlpeps::CompressMPSScheme::VARIATION2Site;
  }
  if (key == "var1" || key == "variational1site" || key == "variation1site" ||
      key == "single-site variational compression" || key == "one_site_variational") {
    return qlpeps::CompressMPSScheme::VARIATION1Site;
  }
  throw std::invalid_argument("MPSCompressScheme must be one of: SVD, Variational2Site, Variational1Site.");
}

/**
 * @brief Common physical parameters for lattice models
 */
struct PhysicalParams : public qlmps::CaseParamsParserBasic {
  PhysicalParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    J2 = ParseDouble("J2");
    RemoveCorner = ParseBool("RemoveCorner");
    // Required: avoid ambiguous defaults in scientific runs.
    // Make users specify explicitly what they simulate.
    try {
      ModelType = ParseStr("ModelType");
    } catch (...) {
      throw std::invalid_argument(
          "Missing required key 'ModelType' in physics params. "
          "Example values: SquareHeisenberg, SquareXY, TriangleHeisenberg.");
    }
    // BoundaryCondition is optional, but must not silently fall back on invalid values.
    // If the key exists and parsing fails, throw to avoid accidental OBC simulations.
    if (this->Has("BoundaryCondition")) {
      BoundaryCondition = ParseBoundaryCondition(ParseStr("BoundaryCondition"));
    } else {
      BoundaryCondition = qlpeps::BoundaryCondition::Open;
    }
  }

  size_t Lx;          ///< Lattice size in x direction
  size_t Ly;          ///< Lattice size in y direction
  double J2;          ///< J1-J2 model parameter (J1 = 1.0 by default)
  bool RemoveCorner;  ///< Whether to remove corner sites (for specific geometries)
  std::string ModelType; ///< Physics model identifier (e.g., SquareHeisenberg, SquareXY, TriangleHeisenberg)
  qlpeps::BoundaryCondition BoundaryCondition; ///< Open or Periodic
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
    MC_total_samples = ParseInt("MC_total_samples");
    WarmUp = ParseInt("WarmUp");
    MCLocalUpdateSweepsBetweenSample = ParseInt("MCLocalUpdateSweepsBetweenSample");
    MCRestrictU1 = ParseBoolOr("MCRestrictU1", true);
  }

  size_t MC_total_samples;                  ///< Total Monte Carlo samples across all MPI ranks
  size_t WarmUp;                           ///< Number of warm-up sweeps
  size_t MCLocalUpdateSweepsBetweenSample; ///< Sweeps between successive samples
  bool MCRestrictU1;                       ///< Whether MC sampler restricts U1 (true by default)
};

/**
 * @brief Boundary MPS specific parameters
 */
struct BMPSParams : public qlmps::CaseParamsParserBasic {
  BMPSParams(const char *f) : CaseParamsParserBasic(f) {
    // Boundary MPS bond dimensions
    // NOTE:
    // These keys are required for OBC(BMPS), but not required for PBC(TRG).
    // We parse them opportunistically here and validate later in the driver based on boundary condition.
    Db_max = this->Has("Dbmps_max") ? static_cast<size_t>(ParseInt("Dbmps_max")) : 0;
    if (this->Has("Dbmps_min")) {
      Db_min = static_cast<size_t>(ParseInt("Dbmps_min"));
    } else {
      Db_min = Db_max;
    }
    if (this->Has("MPSCompressScheme")) {
      MPSCompressScheme = ParseCompressMPSScheme(ParseStr("MPSCompressScheme"));
      has_mps_compress_scheme_ = true;
    } else {
      has_mps_compress_scheme_ = false;
      MPSCompressScheme = qlpeps::CompressMPSScheme::SVD_COMPRESS;
    }
    // Numerical controls specific to BMPS.
    // Default behavior if omitted: trunc_err = 0 (no truncation-by-error).
    TruncErr = ParseDoubleOr("BMPSTruncErr", 0.0);
    ThreadNum = static_cast<size_t>(ParseIntOr("ThreadNum", 1));

    // Variational compression controls (only used when MPSCompressScheme is variational).
    // Defaults keep current behavior (iter_max=10, convergence_tol=TruncErr).
    if (this->Has("BMPSConvergenceTol")) {
      bmps_convergence_tol_ = ParseDouble("BMPSConvergenceTol");
    }
    if (this->Has("BMPSIterMax")) {
      bmps_iter_max_ = static_cast<size_t>(ParseInt("BMPSIterMax"));
    }
  }

  size_t Db_min;                              ///< Minimum boundary MPS bond dimension
  size_t Db_max;                              ///< Maximum boundary MPS bond dimension
  qlpeps::CompressMPSScheme MPSCompressScheme; ///< MPS compression scheme
  double TruncErr;                             ///< Truncation error threshold for boundary MPS
  size_t ThreadNum;                            ///< Number of threads for tensor operations
  bool HasBMPSRequiredKeys() const { return (Db_max > 0); }

  /**
   * @brief Create BMPSTruncatePara from these parameters
   */
  qlpeps::BMPSTruncateParams<qlten::QLTEN_Double> CreateTruncatePara(double trunc_err) const {
    return CreateTruncateParamsImpl_(trunc_err);
  }

  /**
   * @brief Create BMPSTruncatePara using internal TruncErr
   */
  qlpeps::BMPSTruncateParams<qlten::QLTEN_Double> CreateTruncatePara() const {
    return CreateTruncateParamsImpl_(TruncErr);
  }

 private:
  bool has_mps_compress_scheme_ = false;
  std::optional<double> bmps_convergence_tol_;
  std::optional<size_t> bmps_iter_max_;

  qlpeps::BMPSTruncateParams<qlten::QLTEN_Double> CreateTruncateParamsImpl_(double trunc_err) const {
    // Convergence tolerance is an algorithmic knob; don't make it "accidentally zero" when trunc_err=0.
    const double default_tol = (trunc_err > 0.0) ? trunc_err : 1e-12;
    const double tol = bmps_convergence_tol_.value_or(default_tol);
    const size_t it = bmps_iter_max_.value_or(static_cast<size_t>(10));

    switch (MPSCompressScheme) {
      case qlpeps::CompressMPSScheme::SVD_COMPRESS:
        return qlpeps::BMPSTruncateParams<qlten::QLTEN_Double>::SVD(Db_min, Db_max, trunc_err);
      case qlpeps::CompressMPSScheme::VARIATION2Site:
        return qlpeps::BMPSTruncateParams<qlten::QLTEN_Double>::Variational2Site(Db_min, Db_max, trunc_err, tol, it);
      case qlpeps::CompressMPSScheme::VARIATION1Site:
        return qlpeps::BMPSTruncateParams<qlten::QLTEN_Double>::Variational1Site(Db_min, Db_max, trunc_err, tol, it);
      default:
        // Fall back to SVD semantics to avoid undefined behavior if enum extends.
        return qlpeps::BMPSTruncateParams<qlten::QLTEN_Double>::SVD(Db_min, Db_max, trunc_err);
    }
  }
};

/**
 * @brief IO configuration for wavefunction and MC configuration paths
 */
struct IOParams {
  std::string wavefunction_base = "tps";
  std::string configuration_load_dir;
  std::string configuration_dump_dir;

  void Parse(qlmps::CaseParamsParserBasic &parser) {
    wavefunction_base = parser.ParseStrOr("WavefunctionBase", wavefunction_base);
    configuration_load_dir = parser.ParseStrOr("ConfigurationLoadDir", "");
    configuration_dump_dir = parser.ParseStrOr("ConfigurationDumpDir", "");
    if (configuration_load_dir.empty()) configuration_load_dir = wavefunction_base + "final";
    if (configuration_dump_dir.empty()) configuration_dump_dir = wavefunction_base + "final";
  }
};

/**
 * @brief Create PEPSParams selecting TRG (PBC) or BMPS (OBC) backend.
 *
 * For PBC, requires TRGDmin, TRGDmax, TRGTruncErr in the algorithm_parser.
 * For OBC, requires BMPSParams to have valid Dbmps_max.
 */
inline qlpeps::PEPSParams CreatePEPSParams(
    qlpeps::BoundaryCondition bc,
    const BMPSParams &bmps,
    qlmps::CaseParamsParserBasic &algorithm_parser) {
  if (bc == qlpeps::BoundaryCondition::Periodic) {
    if (!(algorithm_parser.Has("TRGDmin") && algorithm_parser.Has("TRGDmax") &&
          algorithm_parser.Has("TRGTruncErr"))) {
      throw std::invalid_argument(
          "PBC requested but TRG params are missing in algorithm JSON. "
          "Require: TRGDmin, TRGDmax, TRGTruncErr (optional: TRGInvRelativeEps).");
    }
    const size_t d_min = static_cast<size_t>(algorithm_parser.ParseInt("TRGDmin"));
    const size_t d_max = static_cast<size_t>(algorithm_parser.ParseInt("TRGDmax"));
    const double trunc_err = algorithm_parser.ParseDouble("TRGTruncErr");
    const double inv_eps = algorithm_parser.ParseDoubleOr("TRGInvRelativeEps", 1e-12);
    return qlpeps::PEPSParams(
        qlpeps::TRGTruncateParams<qlten::QLTEN_Double>(d_min, d_max, trunc_err, inv_eps));
  }

  if (!bmps.HasBMPSRequiredKeys()) {
    throw std::invalid_argument(
        "OBC requested but BMPS params are missing in algorithm JSON. "
        "Require: Dbmps_max (Dbmps_min optional; MPSCompressScheme optional, default=SVD).");
  }
  return qlpeps::PEPSParams(bmps.CreateTruncatePara());
}

/**
 * @brief Parameters for Simple Update (imaginary time evolution)
 */
struct SimpleUpdateParams : public qlmps::CaseParamsParserBasic {
  PhysicalParams physical_params;
  NumericalParams numerical_params;
  double Tau;    ///< Time step for imaginary time evolution
  size_t Step;   ///< Number of simple update steps

  SimpleUpdateParams(const char *physics_file, const char *algorithm_file)
      : qlmps::CaseParamsParserBasic(algorithm_file),
        physical_params(physics_file),
        numerical_params(algorithm_file) {
    Tau = ParseDouble("Tau");
    Step = ParseInt("Step");
  }
};

} // namespace heisenberg_params

#endif //HEISENBERGVMCPEPS_COMMON_PARAMS_H
