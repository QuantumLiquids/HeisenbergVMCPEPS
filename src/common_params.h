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
#include "qlpeps/algorithm/simple_update/simple_update.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "qlpeps/two_dim_tn/common/boundary_condition.h"
#include "qlpeps/two_dim_tn/tensor_network_2d/trg/trg_contractor.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace heisenberg_params {

inline std::string TrimAsciiWhitespace(std::string s) {
  auto is_space = [](unsigned char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
  };
  while (!s.empty() && is_space(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
  while (!s.empty() && is_space(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

inline std::vector<std::string> SplitCsvTokens(
    const std::string &csv,
    const std::string &key_name) {
  const std::string trimmed_csv = TrimAsciiWhitespace(csv);
  if (trimmed_csv.empty()) {
    throw std::invalid_argument(key_name + " must not be empty.");
  }
  if (trimmed_csv.front() == ',' || trimmed_csv.back() == ',') {
    throw std::invalid_argument(key_name + " has malformed leading/trailing comma.");
  }

  std::vector<std::string> tokens;
  std::stringstream stream(trimmed_csv);
  std::string token;
  while (std::getline(stream, token, ',')) {
    token = TrimAsciiWhitespace(token);
    if (token.empty()) {
      throw std::invalid_argument(key_name + " contains an empty item.");
    }
    tokens.push_back(token);
  }
  if (tokens.empty()) {
    throw std::invalid_argument(key_name + " must not be empty.");
  }
  return tokens;
}

inline std::vector<double> ParsePositiveDoubleCsv(
    const std::string &csv,
    const std::string &key_name) {
  const auto tokens = SplitCsvTokens(csv, key_name);
  std::vector<double> values;
  values.reserve(tokens.size());
  for (const auto &token : tokens) {
    size_t parsed_len = 0;
    double value = 0.0;
    try {
      value = std::stod(token, &parsed_len);
    } catch (const std::exception &) {
      throw std::invalid_argument(key_name + " has invalid floating-point item: '" + token + "'.");
    }
    if (parsed_len != token.size()) {
      throw std::invalid_argument(key_name + " has malformed floating-point item: '" + token + "'.");
    }
    if (value <= 0.0) {
      throw std::invalid_argument(key_name + " requires all values > 0.");
    }
    values.push_back(value);
  }
  return values;
}

inline std::vector<size_t> ParsePositiveSizeTCsv(
    const std::string &csv,
    const std::string &key_name) {
  const auto tokens = SplitCsvTokens(csv, key_name);
  std::vector<size_t> values;
  values.reserve(tokens.size());
  for (const auto &token : tokens) {
    size_t parsed_len = 0;
    long long value = 0;
    try {
      value = std::stoll(token, &parsed_len);
    } catch (const std::exception &) {
      throw std::invalid_argument(key_name + " has invalid integer item: '" + token + "'.");
    }
    if (parsed_len != token.size()) {
      throw std::invalid_argument(key_name + " has malformed integer item: '" + token + "'.");
    }
    if (value <= 0) {
      throw std::invalid_argument(key_name + " requires all values > 0.");
    }
    values.push_back(static_cast<size_t>(value));
  }
  return values;
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

enum class InitialConfigStrategy {
  Random,
  Neel,
  ThreeSublatticePolarizedSeed
};

inline InitialConfigStrategy ParseInitialConfigStrategy(const std::string &value) {
  std::string key = TrimAsciiWhitespace(value);
  std::transform(key.begin(), key.end(), key.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (key == "random") return InitialConfigStrategy::Random;
  if (key == "neel") return InitialConfigStrategy::Neel;
  if (key == "threesublatticepolarizedseed") {
    return InitialConfigStrategy::ThreeSublatticePolarizedSeed;
  }
  throw std::invalid_argument(
      "InitialConfigStrategy must be Random, Neel, or ThreeSublatticePolarizedSeed.");
}

/**
 * @brief Common physical parameters for lattice models
 */
struct PhysicalParams : public qlmps::CaseParamsParserBasic {
  PhysicalParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    J2 = ParseDouble("J2");
    RemoveCorner = ParseBoolOr("RemoveCorner", false);
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
  bool RemoveCorner;  ///< Legacy compatibility key. Ignored by unified square/triangle drivers.
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
    initial_config_strategy = ParseInitialConfigStrategy(
        ParseStrOr("InitialConfigStrategy", "Random"));
  }

  size_t MC_total_samples;                  ///< Total Monte Carlo samples across all MPI ranks
  size_t WarmUp;                           ///< Number of warm-up sweeps
  size_t MCLocalUpdateSweepsBetweenSample; ///< Sweeps between successive samples
  bool MCRestrictU1;                       ///< Whether MC sampler restricts U1 (true by default)
  InitialConfigStrategy initial_config_strategy; ///< Init strategy when configuration load fails
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
 * @brief Create initial MC configuration (or load from disk) for half-up/half-down U1 setup.
 *
 * Behavior:
 * - Always starts from a half-up/half-down random configuration.
 * - If load succeeds from configuration_load_dir/configuration{rank}, loaded config wins and warmed_up=true.
 * - If load fails, apply mc_params.initial_config_strategy:
 *   - Random: keep random config.
 *   - Neel: generate checkerboard state with random phase; requires even Lx*Ly.
 *   - ThreeSublatticePolarizedSeed: A sublattice all-up, B/C as close as possible
 *     to quarter-up under exact total Sz=0; requires even Lx*Ly.
 */
inline std::pair<qlpeps::Configuration, bool> InitOrLoadConfigWithStrategy(
    const PhysicalParams &physical_params,
    const MonteCarloNumericalParams &mc_params,
    const std::string &configuration_load_dir,
    int rank) {
  const size_t total_sites = physical_params.Lx * physical_params.Ly;
  const size_t spin_up_sites = total_sites / 2;
  qlpeps::OccupancyNum occupancy(2, 0);
  occupancy[1] = spin_up_sites;
  occupancy[0] = total_sites - spin_up_sites;
  qlpeps::Configuration config(physical_params.Ly, physical_params.Lx, occupancy);

  bool warmed_up = false;
  if (!configuration_load_dir.empty()) {
    warmed_up = config.Load(configuration_load_dir, static_cast<size_t>(rank));
  }
  if (warmed_up) {
    return {config, true};
  }

  if (mc_params.initial_config_strategy == InitialConfigStrategy::Neel) {
    if ((total_sites & 1ULL) != 0ULL) {
      throw std::invalid_argument(
          "InitialConfigStrategy=Neel requires even Lx*Ly when no configuration is loaded.");
    }
    const bool pbc_checkerboard_frustrated =
        (physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic) &&
        (((physical_params.Lx & 1ULL) != 0ULL) || ((physical_params.Ly & 1ULL) != 0ULL));
    if (pbc_checkerboard_frustrated) {
      std::cerr << "[warn][rank " << rank
                << "] InitialConfigStrategy=Neel with PBC and odd Lx or Ly cannot realize a "
                   "true checkerboard AFM on wrap bonds. Proceeding with the checkerboard seed "
                   "pattern anyway."
                << std::endl;
    }
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> phase_dist(0, 1);
    const size_t phase = static_cast<size_t>(phase_dist(rng));
    for (size_t row = 0; row < physical_params.Ly; ++row) {
      for (size_t col = 0; col < physical_params.Lx; ++col) {
        config({row, col}) = (row + col + phase) & 1ULL;
      }
    }
  }
  if (mc_params.initial_config_strategy == InitialConfigStrategy::ThreeSublatticePolarizedSeed) {
    if (physical_params.ModelType != "TriangleHeisenberg") {
      std::cerr << "[warn][rank " << rank
                << "] InitialConfigStrategy=ThreeSublatticePolarizedSeed is used with ModelType="
                << physical_params.ModelType
                << ". Proceeding with the same three-sublattice seed rule." << std::endl;
    }
    if ((total_sites & 1ULL) != 0ULL) {
      throw std::invalid_argument(
          "InitialConfigStrategy=ThreeSublatticePolarizedSeed requires even Lx*Ly when no configuration is loaded.");
    }

    auto sublattice_index = [](size_t row, size_t col) -> size_t {
      return ((row % 3ULL) + 3ULL - (col % 3ULL)) % 3ULL;
    };

    std::vector<std::pair<size_t, size_t>> sub_a_sites;
    std::vector<std::pair<size_t, size_t>> sub_b_sites;
    std::vector<std::pair<size_t, size_t>> sub_c_sites;
    sub_a_sites.reserve(total_sites / 3ULL + 2ULL);
    sub_b_sites.reserve(total_sites / 3ULL + 2ULL);
    sub_c_sites.reserve(total_sites / 3ULL + 2ULL);
    for (size_t row = 0; row < physical_params.Ly; ++row) {
      for (size_t col = 0; col < physical_params.Lx; ++col) {
        const size_t sublattice = sublattice_index(row, col);
        if (sublattice == 0ULL) {
          sub_a_sites.emplace_back(row, col);
        } else if (sublattice == 1ULL) {
          sub_b_sites.emplace_back(row, col);
        } else {
          sub_c_sites.emplace_back(row, col);
        }
      }
    }

    const size_t target_up = total_sites / 2ULL;
    const size_t up_a = sub_a_sites.size();
    const long long remain_ll = static_cast<long long>(target_up) - static_cast<long long>(up_a);
    if (remain_ll < 0 ||
        remain_ll > static_cast<long long>(sub_b_sites.size() + sub_c_sites.size())) {
      std::cerr << "[warn][rank " << rank
                << "] ThreeSublatticePolarizedSeed is infeasible on this lattice; falling back to Random."
                << std::endl;
      return {config, false};
    }

    const size_t remain = static_cast<size_t>(remain_ll);
    const size_t min_up_b = (remain > sub_c_sites.size()) ? (remain - sub_c_sites.size()) : 0ULL;
    const size_t max_up_b = std::min(sub_b_sites.size(), remain);
    if (min_up_b > max_up_b) {
      std::cerr << "[warn][rank " << rank
                << "] ThreeSublatticePolarizedSeed integer assignment failed; falling back to Random."
                << std::endl;
      return {config, false};
    }

    const double ideal_up_b = static_cast<double>(sub_b_sites.size()) * 0.25;
    const double ideal_up_c = static_cast<double>(sub_c_sites.size()) * 0.25;
    double best_cost = std::numeric_limits<double>::infinity();
    size_t best_up_b = min_up_b;
    size_t best_up_c = remain - best_up_b;
    for (size_t up_b = min_up_b; up_b <= max_up_b; ++up_b) {
      const size_t up_c = remain - up_b;
      const double diff_b = static_cast<double>(up_b) - ideal_up_b;
      const double diff_c = static_cast<double>(up_c) - ideal_up_c;
      const double cost = diff_b * diff_b + diff_c * diff_c;
      if (cost < best_cost) {
        best_cost = cost;
        best_up_b = up_b;
        best_up_c = up_c;
      }
    }

    for (size_t row = 0; row < physical_params.Ly; ++row) {
      for (size_t col = 0; col < physical_params.Lx; ++col) {
        config({row, col}) = 0ULL;
      }
    }
    for (const auto &[row, col] : sub_a_sites) {
      config({row, col}) = 1ULL;
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    auto assign_random_up = [&](std::vector<std::pair<size_t, size_t>> &sites, size_t up_count) {
      std::shuffle(sites.begin(), sites.end(), rng);
      for (size_t i = 0; i < up_count; ++i) {
        const auto &[row, col] = sites[i];
        config({row, col}) = 1ULL;
      }
    };
    assign_random_up(sub_b_sites, best_up_b);
    assign_random_up(sub_c_sites, best_up_c);
  }
  return {config, false};
}

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
  /**
   * @brief Optional advanced stop controls for simple update convergence.
   *
   * Semantics:
   * - Convergence gate is (energy criterion) AND (lambda criterion).
   * - Gate must pass for @ref patience consecutive sweeps.
   * - Advanced stop is allowed only after @ref min_steps sweeps.
   */
  struct AdvancedStopOptions {
    double energy_abs_tol = 1e-8;   ///< Absolute energy tolerance.
    double energy_rel_tol = 1e-10;  ///< Relative energy tolerance.
    double lambda_rel_tol = 1e-6;   ///< Relative lambda-drift tolerance.
    size_t patience = 3;            ///< Required consecutive gate passes.
    size_t min_steps = 10;          ///< Minimum executed sweeps before stopping.
  };

  /**
   * @brief Optional multi-stage tau schedule configuration.
   *
   * Semantics:
   * - Stages run in listed order on the same evolving PEPS state.
   * - Each stage uses one tau and one step cap.
   * - If step_caps is omitted in input, all stages inherit global @ref Step.
   */
  struct TauScheduleOptions {
    std::vector<double> taus;       ///< Stage taus in execution order.
    std::vector<size_t> step_caps;  ///< Per-stage step caps.
    bool require_converged = true;  ///< Mark stage failed when summary.converged is false.
    bool dump_each_stage = false;   ///< Whether to dump stage snapshots under dump_dir.
    std::string dump_dir = "tau_schedule";  ///< Driver-owned output directory.
    bool abort_on_stage_failure = true;     ///< Abort schedule on first failed stage.
  };

  PhysicalParams physical_params;
  NumericalParams numerical_params;
  double Tau;    ///< Time step for imaginary time evolution
  size_t Step;   ///< Number of simple update steps
  std::optional<AdvancedStopOptions> advanced_stop;  ///< Advanced stop config; absent means fixed-step mode.
  std::optional<TauScheduleOptions> tau_schedule;    ///< Optional tau schedule; absent means legacy single-stage mode.

  SimpleUpdateParams(const char *physics_file, const char *algorithm_file)
      : qlmps::CaseParamsParserBasic(algorithm_file),
        physical_params(physics_file),
        numerical_params(algorithm_file) {
    Tau = ParseDouble("Tau");
    const int parsed_step = ParseInt("Step");
    Step = static_cast<size_t>(parsed_step);

    const bool has_advanced_stop_enabled_key = this->Has("AdvancedStopEnabled");
    const bool advanced_stop_enabled = ParseBoolOr("AdvancedStopEnabled", false);
    const bool has_advanced_stop_tuning_key =
        this->Has("AdvancedStopEnergyAbsTol") ||
        this->Has("AdvancedStopEnergyRelTol") ||
        this->Has("AdvancedStopLambdaRelTol") ||
        this->Has("AdvancedStopPatience") ||
        this->Has("AdvancedStopMinSteps");

    bool enable_advanced_stop = false;
    if (has_advanced_stop_enabled_key && !advanced_stop_enabled) {
      // Explicit false wins over implicit auto-enable from tuning keys.
      enable_advanced_stop = false;
    } else {
      enable_advanced_stop = advanced_stop_enabled || has_advanced_stop_tuning_key;
    }

    if (enable_advanced_stop) {
      AdvancedStopOptions options;
      options.energy_abs_tol = ParseDoubleOr("AdvancedStopEnergyAbsTol", options.energy_abs_tol);
      options.energy_rel_tol = ParseDoubleOr("AdvancedStopEnergyRelTol", options.energy_rel_tol);
      options.lambda_rel_tol = ParseDoubleOr("AdvancedStopLambdaRelTol", options.lambda_rel_tol);

      const int patience = ParseIntOr("AdvancedStopPatience", static_cast<int>(options.patience));
      const int min_steps = ParseIntOr("AdvancedStopMinSteps", static_cast<int>(options.min_steps));

      if (options.energy_abs_tol < 0.0) {
        throw std::invalid_argument("AdvancedStopEnergyAbsTol must be >= 0.");
      }
      if (options.energy_rel_tol < 0.0) {
        throw std::invalid_argument("AdvancedStopEnergyRelTol must be >= 0.");
      }
      if (options.lambda_rel_tol < 0.0) {
        throw std::invalid_argument("AdvancedStopLambdaRelTol must be >= 0.");
      }
      if (patience <= 0) {
        throw std::invalid_argument("AdvancedStopPatience must be > 0.");
      }
      if (min_steps <= 0) {
        throw std::invalid_argument("AdvancedStopMinSteps must be > 0.");
      }

      options.patience = static_cast<size_t>(patience);
      options.min_steps = static_cast<size_t>(min_steps);
      advanced_stop = options;
    }

    const bool tau_schedule_enabled = ParseBoolOr("TauScheduleEnabled", false);
    if (tau_schedule_enabled) {
      TauScheduleOptions options;
      if (!this->Has("TauScheduleTaus")) {
        throw std::invalid_argument(
            "TauScheduleEnabled=true requires key TauScheduleTaus.");
      }
      options.taus = ParsePositiveDoubleCsv(ParseStr("TauScheduleTaus"), "TauScheduleTaus");
      if (options.taus.empty()) {
        throw std::invalid_argument(
            "TauScheduleEnabled=true requires non-empty TauScheduleTaus.");
      }

      if (this->Has("TauScheduleStepCaps")) {
        options.step_caps = ParsePositiveSizeTCsv(
            ParseStr("TauScheduleStepCaps"),
            "TauScheduleStepCaps");
        if (options.step_caps.size() != options.taus.size()) {
          throw std::invalid_argument(
              "TauScheduleStepCaps count must match TauScheduleTaus count.");
        }
      } else {
        if (parsed_step <= 0) {
          throw std::invalid_argument(
              "Step must be > 0 when TauScheduleEnabled=true and TauScheduleStepCaps is omitted.");
        }
        options.step_caps.assign(options.taus.size(), static_cast<size_t>(parsed_step));
      }

      options.require_converged = ParseBoolOr("TauScheduleRequireConverged", true);
      if (options.require_converged && !advanced_stop.has_value()) {
        throw std::invalid_argument(
            "TauScheduleRequireConverged=true requires advanced stop to be enabled. "
            "Set AdvancedStopEnabled=true (or provide AdvancedStop* keys), "
            "or set TauScheduleRequireConverged=false.");
      }
      options.dump_each_stage = ParseBoolOr("TauScheduleDumpEachStage", false);
      options.dump_dir = TrimAsciiWhitespace(ParseStrOr("TauScheduleDumpDir", "tau_schedule"));
      if (options.dump_dir.empty()) {
        throw std::invalid_argument("TauScheduleDumpDir must not be empty.");
      }
      options.abort_on_stage_failure = ParseBoolOr("TauScheduleAbortOnStageFailure", true);
      tau_schedule = options;
    }
  }

  /**
   * @brief Build qlpeps simple-update parameters from parsed driver inputs.
   */
  qlpeps::SimpleUpdatePara CreateSimpleUpdatePara() const {
    return CreateSimpleUpdateParaForStage(Tau, Step);
  }

  /**
   * @brief Build per-stage qlpeps simple-update parameters for tau-schedule mode.
   */
  qlpeps::SimpleUpdatePara CreateSimpleUpdateParaForStage(
      double tau,
      size_t step_cap) const {
    if (!advanced_stop.has_value()) {
      return qlpeps::SimpleUpdatePara(
          step_cap,
          tau,
          numerical_params.Dmin,
          numerical_params.Dmax,
          numerical_params.TruncErr);
    }
    const auto &cfg = advanced_stop.value();
    return qlpeps::SimpleUpdatePara::Advanced(
        step_cap,
        tau,
        numerical_params.Dmin,
        numerical_params.Dmax,
        numerical_params.TruncErr,
        cfg.energy_abs_tol,
        cfg.energy_rel_tol,
        cfg.lambda_rel_tol,
        cfg.patience,
        cfg.min_steps);
  }
};

} // namespace heisenberg_params

#endif //HEISENBERGVMCPEPS_COMMON_PARAMS_H
