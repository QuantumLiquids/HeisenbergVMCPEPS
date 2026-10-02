// SPDX-License-Identifier: LGPL-3.0-only
//
// Runtime check of the opt-in axis update on a state of this repository. The updater that
// VisitOBCUpdater selects for MCAxisUpdate=true,
// MCUpdateSequence<MCUpdateSquareAxisAutoregressiveOBC, MCUpdateSquareTNN3SiteExchangeOBC>,
// runs Monte Carlo sweeps on an Lx = 4, Ly = 3 Heisenberg state of this build's types
// (qldouble.h: TenElemT, QNT, pb_out, label 0 = up), made the way simple_update makes one (Neel
// start, Heisenberg bond gate, TPS tensors normalized). The parameters come from algorithm JSON
// files through the VMC and measurement parsers, and the BMPS truncation through the drivers'
// CreateContractorParams. The sweeps run directly on a TPSWaveFunctionComponent: no MPI and no
// VMC or measurement run.
//
// Checks: N_up is conserved by every sweep when S_z is conserved (dense tensors with
// MCRestrictU1=true, which maps to the N_up count table, or -DU1SYM tensors); slices move; the
// component amplitude equals a cold contraction (the BMPS is exact here); no zero-support visit;
// the bound floor(max(Lx, Ly) / 2) + 1 that the parsers enforce is reached on a real state and is
// tight (one less makes the library throw at the first row visit); the accept_rates and
// [MC updater] entries the tutorial documents.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>
#include "qlpeps/algorithm/simple_update/square_lattice_nn_simple_update.h"
#include "qlpeps/api/conversions.h"
#include "../src/enhanced_params_parser.h"
#include "../src/enhanced_measure_params_parser.h"
#include "../src/model_updater_factory.h"
#include "../src/qldouble.h"

namespace axis_update_sweep_test {

using qlpeps::config::Json;
using LocalUpdater = qlpeps::MCUpdateSquareTNN3SiteExchangeOBC;
using AxisUpdater = qlpeps::MCUpdateSquareAxisAutoregressiveOBC;
using Sequence = qlpeps::MCUpdateSequence<AxisUpdater, LocalUpdater>;
using Component = qlpeps::TPSWaveFunctionComponent<TenElemT, QNT>;
using SITPS = qlpeps::SplitIndexTPS<TenElemT, QNT>;

constexpr bool kDense = std::is_same_v<QNT, qlten::special_qn::TrivialRepQN>;

// Lattice: Ly = 3 rows of length Lx = 4, so rows are the longer slices and set the bound
// floor(max(Lx, Ly) / 2) + 1 = 3; columns (length 3) have at most 2 sector blocks per bond. Under
// the Neel configuration every row has N_up = 2, so the middle bond of row 0, the first slice a
// sweep visits, carries the 3 sector blocks N_up = 0, 1, 2 of its first two sites before any
// draw has changed the configuration.
constexpr size_t kLx = 4;
constexpr size_t kLy = 3;
constexpr size_t kBound = 3;
constexpr size_t kPEPSBondDim = 4;
constexpr size_t kBMPSBondDim = 16;  // >= 4^2: the boundary MPS of this lattice is exact
constexpr size_t kSweeps = 4;
constexpr std::uint32_t kSeed = 20261002;

void Check(bool ok, const char *expression, int line) {
  if (!ok) {
    throw std::runtime_error("test_axis_update_sweep.cpp:" + std::to_string(line) +
                             ": check failed: " + expression);
  }
}
#define CHECK(condition) Check((condition), #condition, __LINE__)

bool Contains(const std::string &text, const std::string &part) {
  return text.find(part) != std::string::npos;
}

/// Message of the std::invalid_argument @p f throws; empty when it throws nothing.
template<typename F>
std::string InvalidArgumentMessage(F &&f) {
  try {
    f();
  } catch (const std::invalid_argument &error) {
    return error.what();
  }
  return "";
}

/// Temporary physics/algorithm JSON files, removed with their directory.
class TempFiles {
 public:
  TempFiles()
      : dir_(std::filesystem::temp_directory_path() /
             ("axis-update-sweep-" +
              std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))) {
    std::filesystem::create_directories(dir_);
  }
  ~TempFiles() {
    std::error_code ignored;
    std::filesystem::remove_all(dir_, ignored);
  }
  TempFiles(const TempFiles &) = delete;
  TempFiles &operator=(const TempFiles &) = delete;

  /// Write {"CaseParams": case_params} to @p name and return its path.
  std::string Write(const std::string &name, const Json &case_params) const {
    const auto path = dir_ / name;
    std::ofstream out(path);
    out << Json{{"CaseParams", case_params}};
    out.close();
    return path.string();
  }

 private:
  std::filesystem::path dir_;
};

/// Neel labels as simple_update seeds them: label (row + col) % 2.
qlpeps::Configuration NeelConfiguration() {
  qlpeps::Configuration config(kLy, kLx);
  for (size_t row = 0; row < kLy; ++row) {
    for (size_t col = 0; col < kLx; ++col) {
      config({row, col}) = (row + col) % 2;
    }
  }
  return config;
}

size_t CountUp(const qlpeps::Configuration &config) {
  size_t n_up = 0;
  for (size_t row = 0; row < kLy; ++row) {
    for (size_t col = 0; col < kLx; ++col) {
      n_up += (config({row, col}) == 0);  // label 0 = up (qldouble.h)
    }
  }
  return n_up;
}

/**
 * The Heisenberg state of this build's types, made as simple_update makes it: Neel product PEPS,
 * the same nearest-neighbour gate, 10 steps at tau = 0.1 (D = 4), then ToTPS, every tensor
 * normalized, and ToSplitIndexTPS. Its Neel amplitude is nonzero by construction.
 */
SITPS HeisenbergState() {
  Tensor ham_hei_nn = Tensor({pb_in, pb_out, pb_in, pb_out});
  ham_hei_nn({0, 0, 0, 0}) = 0.25;
  ham_hei_nn({1, 1, 1, 1}) = 0.25;
  ham_hei_nn({1, 1, 0, 0}) = -0.25;
  ham_hei_nn({0, 0, 1, 1}) = -0.25;
  ham_hei_nn({0, 1, 1, 0}) = 0.5;
  ham_hei_nn({1, 0, 0, 1}) = 0.5;

  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps0(pb_out, kLy, kLx, qlpeps::BoundaryCondition::Open);
  std::vector<std::vector<size_t>> activates(kLy, std::vector<size_t>(kLx));
  for (size_t row = 0; row < kLy; ++row) {
    for (size_t col = 0; col < kLx; ++col) {
      activates[row][col] = (row + col) % 2;
    }
  }
  peps0.Initial(activates);

  const qlpeps::SimpleUpdateParams update_params(/*steps=*/10, /*tau=*/0.1, kPEPSBondDim,
                                                 kPEPSBondDim, /*trunc_err=*/1e-10);
  qlpeps::SquareLatticePEPS<TenElemT, QNT> peps(pb_out, kLy, kLx, qlpeps::BoundaryCondition::Open);
  {
    // The executor prints a banner when constructed and a line per step; discard both.
    std::ostringstream discarded;
    struct RestoreCout {
      std::streambuf *previous;
      ~RestoreCout() { std::cout.rdbuf(previous); }
    } restore{std::cout.rdbuf(discarded.rdbuf())};
    qlpeps::SquareLatticeNNSimpleUpdateExecutor<TenElemT, QNT> executor(update_params, peps0,
                                                                         ham_hei_nn);
    executor.Execute();
    peps = executor.GetPEPS();
  }

  auto tps = qlpeps::ToTPS<TenElemT, QNT>(peps);
  for (auto &tensor : tps) {
    tensor.Normalize();
  }
  return qlpeps::ToSplitIndexTPS<TenElemT, QNT>(tps);
}

Json Physics() {
  return Json{{"Lx", kLx}, {"Ly", kLy}, {"J2", 0.0}, {"ModelType", "SquareHeisenberg"},
              {"BoundaryCondition", "Open"}};
}

/// An OBC algorithm file with the axis update on and an exact BMPS.
Json Algorithm() {
  return Json{{"MC_total_samples", 8}, {"WarmUp", 0}, {"MCLocalUpdateSweepsBetweenSample", 1},
              {"Dbmps_max", kBMPSBondDim}, {"OptimizerType", "SGD"}, {"MCAxisUpdate", true}};
}

/// The drivers' BMPS truncation (heisenberg_params::CreateContractorParams, OBC branch).
qlpeps::BMPSTruncateParams<double> TruncateParams(const heisenberg_params::BMPSParams &bmps) {
  return heisenberg_params::CreateContractorParams(qlpeps::BoundaryCondition::Open, bmps,
                                                   bmps.algorithm_values)
      .Get<qlpeps::ContractorParams::BMPSParams>();
}

/// What a run of sweeps showed.
struct SweepResult {
  qlpeps::AxisSweepStats axis;  ///< cumulative over the run
  std::vector<std::pair<std::string, double>> diagnostics;  ///< the sequence's Diagnostics()
  bool count_table = false;     ///< the axis child has a count table
};

/**
 * Run kSweeps sweeps of the updater VisitOBCUpdater selects, from the Neel configuration, and
 * check after every sweep: the updater is the sequence; the accept_rates are {changed rows / Ly,
 * changed columns / Lx, local acceptance}; N_up is unchanged when @p conserves_sz; and the
 * component amplitude equals a cold contraction of its configuration.
 */
SweepResult RunSweeps(const SITPS &sitps,
                      const heisenberg_params::AxisUpdateParams &axis,
                      const heisenberg_params::MonteCarloNumericalParams &mc,
                      const qlpeps::BMPSTruncateParams<double> &trunc,
                      bool conserves_sz) {
  SweepResult result;
  size_t visits = 0;
  VisitOBCUpdater<QNT>(axis, mc, [&](auto updater) {
    ++visits;
    if constexpr (!std::is_same_v<decltype(updater), Sequence>) {
      throw std::runtime_error("VisitOBCUpdater did not select the axis sequence");
    } else {
      updater.SeedRandomEngine(kSeed);
      result.count_table = updater.template Get<0>().Params().count_constraint.has_value();
      Component component(sitps, NeelConfiguration(), trunc);
      const size_t n_up = CountUp(component.config);
      CHECK(n_up == kLx * kLy / 2);
      for (size_t sweep = 0; sweep < kSweeps; ++sweep) {
        std::vector<double> accept_rates;
        updater(sitps, component, accept_rates);
        const qlpeps::AxisSweepStats &last = updater.template Get<0>().LastSweepStats();
        CHECK(last.row_visits == kLy && last.col_visits == kLx);
        CHECK(accept_rates.size() == 3);
        CHECK(accept_rates[0] == double(last.rows_changed) / double(kLy));
        CHECK(accept_rates[1] == double(last.cols_changed) / double(kLx));
        CHECK(accept_rates[2] >= 0.0 && accept_rates[2] <= 1.0);
        if (conserves_sz) {
          CHECK(CountUp(component.config) == n_up);
        }
        const Component cold(sitps, component.config, trunc);
        CHECK(std::isfinite(std::abs(cold.amplitude)) && std::abs(cold.amplitude) > 0.0);
        CHECK(std::abs(component.amplitude - cold.amplitude) <= 1e-8 * std::abs(cold.amplitude));
      }
      result.axis = updater.template Get<0>().CumulativeStats();
      for (const auto &entry : updater.Diagnostics()) {
        result.diagnostics.emplace_back(entry.name, entry.value);
      }
    }
  });
  CHECK(visits == 1);
  return result;
}

double DiagnosticValue(const SweepResult &result, const std::string &name) {
  for (const auto &[entry, value] : result.diagnostics) {
    if (entry == name) return value;
  }
  throw std::runtime_error("no diagnostic named " + name);
}

/// Statistics and [MC updater] entries every run must show.
void CheckRunStatistics(const SweepResult &result, size_t expected_block_count) {
  CHECK(result.axis.row_visits == kSweeps * kLy);
  CHECK(result.axis.col_visits == kSweeps * kLx);
  CHECK(result.axis.rows_changed + result.axis.cols_changed > 0);  // the slices move
  CHECK(result.axis.zero_support_old == 0);
  CHECK(result.axis.max_block_count == expected_block_count);

  // The [MC updater] entries the tutorial names: child0.axis.*, only from the axis child.
  const std::set<std::string> expected = {
      "child0.axis.row_visits", "child0.axis.rows_changed", "child0.axis.col_visits",
      "child0.axis.cols_changed", "child0.axis.zero_support_old",
      "child0.axis.max_discarded_weight", "child0.axis.max_slice_bond_dim",
      "child0.axis.max_block_count"};
  std::set<std::string> names;
  for (const auto &entry : result.diagnostics) names.insert(entry.first);
  CHECK(names == expected);
  CHECK(DiagnosticValue(result, "child0.axis.row_visits") == double(kSweeps * kLy));
  CHECK(DiagnosticValue(result, "child0.axis.zero_support_old") == 0.0);
  CHECK(DiagnosticValue(result, "child0.axis.max_block_count") == double(expected_block_count));
}

/**
 * S_z conserved: dense tensors with MCRestrictU1=true (the N_up table), or U1QN tensors with either
 * MCRestrictU1 value (no table). MCAxisUpdateDmax at the parsers' bound; the measurement parser
 * reads the file when @p measure_parser.
 */
void TestConservingSweeps(const SITPS &sitps, const TempFiles &files, bool mc_restrict_u1,
                          bool measure_parser) {
  const auto physics = files.Write("physics.json", Physics());
  Json values = Algorithm();
  values["MCRestrictU1"] = mc_restrict_u1;
  values["MCAxisUpdateDmax"] = kBound;
  const auto algorithm = files.Write("conserving.json", values);
  const auto run = [&](const auto &params) {
    CHECK(heisenberg_params::AxisUpdateParams::SectorBlockBound(params.physical_params) == kBound);
    const heisenberg_params::AxisUpdateParams &axis = params.axis_update_params;
    CHECK(axis.MCAxisUpdate && axis.MCAxisUpdateDmax == kBound);
    return RunSweeps(sitps, axis, params.mc_params, TruncateParams(params.bmps_params),
                     /*conserves_sz=*/true);
  };
  const SweepResult result = measure_parser
      ? run(EnhancedMCMeasureParams(physics.c_str(), algorithm.c_str()))
      : run(EnhancedVMCUpdateParams(physics.c_str(), algorithm.c_str()));
  CHECK(result.count_table == kDense);  // the table exactly for dense tensors
  // The middle bond of a row reaches the bound: the bound is attained, and enough.
  CheckRunStatistics(result, kBound);
}

/**
 * One below the bound (bypassing the parser, which rejects it): the library throws at the first
 * visit, row 0 of the Neel configuration, whose middle bond has 3 sector blocks. So the parsers'
 * bound is tight on this state; the shorter columns alone would allow 2.
 */
void TestBelowTheBoundThrowsInTheLibrary(const SITPS &sitps, const TempFiles &files) {
  const auto physics = files.Write("physics.json", Physics());
  Json values = Algorithm();
  values["MCAxisUpdateDmax"] = kBound;
  const auto algorithm = files.Write("below.json", values);
  const EnhancedVMCUpdateParams params(physics.c_str(), algorithm.c_str());
  heisenberg_params::AxisUpdateParams below = params.axis_update_params;
  below.MCAxisUpdateDmax = kBound - 1;
  CHECK(Contains(InvalidArgumentMessage([&] {
                   below.Validate(params.physical_params, params.mc_params.MCRestrictU1,
                                  /*spin_inversion_projected=*/false,
                                  heisenberg_params::kU1SymmetricBuild);
                 }),
                 "MCAxisUpdateDmax = 2 is below floor(max(Lx, Ly) / 2) + 1 = 3"));
  const std::string message = InvalidArgumentMessage([&] {
    RunSweeps(sitps, below, params.mc_params, TruncateParams(params.bmps_params),
              /*conserves_sz=*/true);
  });
  CHECK(Contains(message, "HORIZONTAL slice 0"));
  CHECK(Contains(message, "has 3 nonzero sector blocks but the slice D_max is 2"));
}

/**
 * Dense tensors with MCRestrictU1=false: no table, so one sector block per bond and no S_z
 * bound; MCAxisUpdateDmax defaults to Dbmps_max (untruncated slice MPS here). N_up is not
 * asserted: the axis moves may change it (tutorial section 4.8).
 */
void TestDenseWithoutTableSweeps(const SITPS &sitps, const TempFiles &files) {
  const auto physics = files.Write("physics.json", Physics());
  Json values = Algorithm();
  values["MCRestrictU1"] = false;
  const auto algorithm = files.Write("free.json", values);
  const EnhancedVMCUpdateParams params(physics.c_str(), algorithm.c_str());
  CHECK(params.axis_update_params.MCAxisUpdateDmax == kBMPSBondDim);
  const SweepResult result = RunSweeps(sitps, params.axis_update_params, params.mc_params,
                                       TruncateParams(params.bmps_params), /*conserves_sz=*/false);
  CHECK(!result.count_table);
  CheckRunStatistics(result, 1);
  CHECK(result.axis.max_discarded_weight == 0.0);  // untruncated slice MPS
}

}  // namespace axis_update_sweep_test

int main() {
  using namespace axis_update_sweep_test;
  // One numerics thread, as the PEPS unit tests pin it: tiny tensors, reproducible rounding.
  qlten::hp_numeric::SetCpuNumericsThreads(1);
  try {
    const TempFiles files;
    const SITPS sitps = HeisenbergState();
    CHECK(sitps.GetBoundaryCondition() == qlpeps::BoundaryCondition::Open);
    CHECK(sitps({0, 0}).size() == 2);  // two labels per site, both covered by the N_up table
    TestConservingSweeps(sitps, files, /*mc_restrict_u1=*/true, /*measure_parser=*/false);
    TestConservingSweeps(sitps, files, /*mc_restrict_u1=*/true, /*measure_parser=*/true);
    if constexpr (kDense) {
      TestDenseWithoutTableSweeps(sitps, files);
    } else {
      TestConservingSweeps(sitps, files, /*mc_restrict_u1=*/false, /*measure_parser=*/false);
    }
    TestBelowTheBoundThrowsInTheLibrary(sitps, files);
  } catch (const std::exception &error) {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  std::cout << "test_axis_update_sweep: all checks passed" << std::endl;
  return 0;
}
