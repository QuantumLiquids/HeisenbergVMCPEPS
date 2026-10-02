// SPDX-License-Identifier: LGPL-3.0-only
//
// Regression for the opt-in axis update keys (MCAxisUpdate, MCAxisUpdateDmin, MCAxisUpdateDmax,
// MCAxisUpdateTruncErr): parsing in the VMC and measurement parsers, validation, the PEPS
// parameters they map to, the updater VisitOBCUpdater selects and the sampler log.
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include "../src/enhanced_params_parser.h"
#include "../src/enhanced_measure_params_parser.h"
#include "../src/model_updater_factory.h"
#include "../src/qldouble.h"

// The parsers' build flag and the drivers' QNT come from the same -DU1SYM switch.
static_assert(heisenberg_params::kU1SymmetricBuild == std::is_same_v<QNT, qlten::special_qn::U1QN>,
              "kU1SymmetricBuild must match the QNT of qldouble.h");

namespace axis_update_params_test {

using qlpeps::config::Json;
using Dense = qlten::special_qn::TrivialRepQN;
using U1 = qlten::special_qn::U1QN;
using LocalUpdater = qlpeps::MCUpdateSquareTNN3SiteExchangeOBC;
using AxisUpdater = qlpeps::MCUpdateSquareAxisAutoregressiveOBC;
using Sequence = qlpeps::MCUpdateSequence<AxisUpdater, LocalUpdater>;

/// The MCRestrictU1=false notice exactly as the drivers printed it before the axis update existed.
const char *const kLegacyRestrictU1Notice =
    "[info] MCRestrictU1=false (no-U1 sampler requested). Using U1 sampler temporarily; "
    "non-U1 variant will be wired in a later step.\n";

void Check(bool ok, const char *expression, int line) {
  if (!ok) {
    throw std::runtime_error("test_axis_update_params.cpp:" + std::to_string(line) +
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

/// Everything @p f writes to std::cout.
template<typename F>
std::string CaptureStdout(F &&f) {
  std::ostringstream buffer;
  std::streambuf *const previous = std::cout.rdbuf(buffer.rdbuf());
  try {
    f();
  } catch (...) {
    std::cout.rdbuf(previous);
    throw;
  }
  std::cout.rdbuf(previous);
  return buffer.str();
}

/// Temporary physics/algorithm JSON files, removed with their directory.
class TempFiles {
 public:
  TempFiles()
      : dir_(std::filesystem::temp_directory_path() /
             ("axis-update-params-" +
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

Json Physics(size_t lx, size_t ly, const std::string &boundary = "Open") {
  return Json{{"Lx", lx}, {"Ly", ly}, {"J2", 0.0}, {"ModelType", "SquareHeisenberg"},
              {"BoundaryCondition", boundary}};
}

/// A minimal OBC algorithm file without axis keys (Dbmps_max = 10).
Json Algorithm() {
  return Json{{"MC_total_samples", 8}, {"WarmUp", 0}, {"MCLocalUpdateSweepsBetweenSample", 1},
              {"Dbmps_max", 10}, {"OptimizerType", "SGD"}};
}

/// Which updater VisitOBCUpdater passed, and the axis child's parameters for a sequence.
struct Visited {
  size_t calls = 0;
  bool local = false;
  bool sequence = false;
  std::optional<qlpeps::AxisAutoregressiveUpdateParams> axis;
};

template<typename QNT>
Visited Visit(const heisenberg_params::AxisUpdateParams &axis,
              const heisenberg_params::MonteCarloNumericalParams &mc) {
  Visited visited;
  VisitOBCUpdater<QNT>(axis, mc, [&](auto updater) {
    using U = decltype(updater);
    ++visited.calls;
    if constexpr (std::is_same_v<U, LocalUpdater>) {
      visited.local = true;
    } else if constexpr (std::is_same_v<U, Sequence>) {
      visited.sequence = true;
      visited.axis = updater.template Get<0>().Params();
    } else {
      static_assert(sizeof(U) == 0, "VisitOBCUpdater passed an unexpected updater type");
    }
  });
  return visited;
}

void CheckSliceParams(const qlpeps::AxisAutoregressiveUpdateParams &params,
                      size_t d_min, size_t d_max, double trunc_err) {
  CHECK(params.slice_params == qlpeps::BMPSSliceSamplerParams(d_min, d_max, trunc_err));
  CHECK(params.axes == qlpeps::SliceAxes::kRowsAndColumns);
  CHECK(!params.exactness_self_check);
}

const qlpeps::LabelCountTable kNUpTable({{1}, {0}});

void TestDefaultsKeepTheLocalUpdater(const TempFiles &files) {
  const auto physics = files.Write("physics.json", Physics(6, 4));
  const auto algorithm = files.Write("algorithm.json", Algorithm());
  const EnhancedVMCUpdateParams vmc(physics.c_str(), algorithm.c_str());
  const EnhancedMCMeasureParams measure(physics.c_str(), algorithm.c_str());
  for (const auto *axis : {&vmc.axis_update_params, &measure.axis_update_params}) {
    CHECK(!axis->MCAxisUpdate);
    CHECK(axis->MCAxisUpdateDmin == 1);
    CHECK(axis->MCAxisUpdateDmax == 10);  // Dbmps_max
    CHECK(axis->dmax_from_dbmps_max);
    CHECK(axis->MCAxisUpdateTruncErr == 0.0);
  }
  for (const Visited &visited : {Visit<Dense>(vmc.axis_update_params, vmc.mc_params),
                                 Visit<U1>(vmc.axis_update_params, vmc.mc_params)}) {
    CHECK(visited.calls == 1);
    CHECK(visited.local);
    CHECK(!visited.sequence);
  }
  // A default-constructed AxisUpdateParams is the disabled state too.
  CHECK(Visit<Dense>(heisenberg_params::AxisUpdateParams{}, vmc.mc_params).local);

  // Logs: nothing with MCRestrictU1=true, on every rank.
  for (int rank : {0, 1}) {
    CHECK(CaptureStdout([&] { LogSamplerChoice<Dense>(vmc.mc_params, vmc.axis_update_params, rank); })
              .empty());
  }

  // Logs: the historical notice, byte for byte, on every rank with MCRestrictU1=false.
  Json no_restrict = Algorithm();
  no_restrict["MCRestrictU1"] = false;
  no_restrict["MCAxisUpdate"] = false;  // explicit false behaves like an absent key
  const auto no_restrict_path = files.Write("no_restrict.json", no_restrict);
  const EnhancedVMCUpdateParams vmc_no_restrict(physics.c_str(), no_restrict_path.c_str());
  CHECK(!vmc_no_restrict.axis_update_params.MCAxisUpdate);
  CHECK(Visit<Dense>(vmc_no_restrict.axis_update_params, vmc_no_restrict.mc_params).local);
  for (int rank : {0, 1}) {
    CHECK(CaptureStdout([&] {
            LogSamplerChoice<Dense>(vmc_no_restrict.mc_params, vmc_no_restrict.axis_update_params, rank);
          }) == kLegacyRestrictU1Notice);
  }
}

void TestExplicitValuesAndMapping(const TempFiles &files) {
  const auto physics = files.Write("physics.json", Physics(6, 4));
  Json values = Algorithm();
  values["MCAxisUpdate"] = true;
  values["MCAxisUpdateDmin"] = 2;
  values["MCAxisUpdateDmax"] = 9;
  values["MCAxisUpdateTruncErr"] = 1e-7;
  const auto algorithm = files.Write("algorithm.json", values);
  const EnhancedVMCUpdateParams vmc(physics.c_str(), algorithm.c_str());
  const EnhancedMCMeasureParams measure(physics.c_str(), algorithm.c_str());
  for (const auto *axis : {&vmc.axis_update_params, &measure.axis_update_params}) {
    CHECK(axis->MCAxisUpdate);
    CHECK(axis->MCAxisUpdateDmin == 2);
    CHECK(axis->MCAxisUpdateDmax == 9);
    CHECK(!axis->dmax_from_dbmps_max);
    CHECK(axis->MCAxisUpdateTruncErr == 1e-7);
  }

  // Dense tensors with MCRestrictU1=true (the default): the sequence, with the N_up table.
  const Visited dense = Visit<Dense>(vmc.axis_update_params, vmc.mc_params);
  CHECK(dense.calls == 1);
  CHECK(dense.sequence && !dense.local);
  CHECK(dense.axis.has_value());
  CheckSliceParams(*dense.axis, 2, 9, 1e-7);
  CHECK(dense.axis->count_constraint.has_value());
  CHECK(*dense.axis->count_constraint == kNUpTable);

  // U1QN tensors carry S_z: no table.
  const Visited u1 = Visit<U1>(vmc.axis_update_params, vmc.mc_params);
  CHECK(u1.sequence);
  CheckSliceParams(*u1.axis, 2, 9, 1e-7);
  CHECK(!u1.axis->count_constraint.has_value());

  // The direct mapping agrees with the visited updater.
  CHECK(MakeAxisAutoregressiveUpdateParams<Dense>(vmc.axis_update_params, true).count_constraint ==
        std::optional<qlpeps::LabelCountTable>(kNUpTable));
  CHECK(!MakeAxisAutoregressiveUpdateParams<Dense>(vmc.axis_update_params, false).count_constraint);
  CHECK(!MakeAxisAutoregressiveUpdateParams<U1>(vmc.axis_update_params, true).count_constraint);
  CHECK(!MakeAxisAutoregressiveUpdateParams<U1>(vmc.axis_update_params, false).count_constraint);

  // Logs: rank 0 only, with the slice MPS and the table choice.
  const std::string dense_log =
      CaptureStdout([&] { LogSamplerChoice<Dense>(vmc.mc_params, vmc.axis_update_params, 0); });
  CHECK(Contains(dense_log, "MCAxisUpdate=true"));
  CHECK(Contains(dense_log, "D_min=2 D_max=9 trunc_err=1e-07"));
  CHECK(Contains(dense_log, "count table: N_up {{1},{0}}"));
  CHECK(Contains(dense_log, "0 = changed rows / Ly, 1 = changed columns / Lx, 2 = local"));
  CHECK(!Contains(dense_log, "no-U1 sampler"));
  CHECK(CaptureStdout([&] { LogSamplerChoice<Dense>(vmc.mc_params, vmc.axis_update_params, 1); })
            .empty());
  CHECK(Contains(CaptureStdout([&] { LogSamplerChoice<U1>(vmc.mc_params, vmc.axis_update_params, 0); }),
                 "count table: none (U1QN tensors conserve S_z)"));

  // Dense tensors with MCRestrictU1=false: no table, so the axis moves change S_z. Without the
  // table there is one sector block per bond, so D_max = 1 is accepted for dense tensors.
  values["MCRestrictU1"] = false;
  values["MCAxisUpdateDmin"] = 1;
  values["MCAxisUpdateDmax"] = 1;
  const auto no_restrict = files.Write("no_restrict.json", values);
  const heisenberg_params::PhysicalParams physical(physics.c_str());
  const heisenberg_params::MonteCarloNumericalParams mc(no_restrict.c_str());
  const heisenberg_params::BMPSParams bmps(no_restrict.c_str());
  const heisenberg_params::AxisUpdateParams dense_free(no_restrict.c_str(), physical, mc, bmps,
                                                       /*spin_inversion_projected=*/false,
                                                       /*u1_symmetric_build=*/false);
  const Visited free_sector = Visit<Dense>(dense_free, mc);
  CHECK(free_sector.sequence);
  CheckSliceParams(*free_sector.axis, 1, 1, 1e-7);
  CHECK(!free_sector.axis->count_constraint.has_value());
  const std::string free_log = CaptureStdout([&] { LogSamplerChoice<Dense>(mc, dense_free, 0); });
  CHECK(Contains(free_log, "MCRestrictU1=false: axis moves change S_z"));
  CHECK(!Contains(free_log, "no-U1 sampler"));
}

void TestDmaxDefaultsToDbmpsMax(const TempFiles &files) {
  const auto physics = files.Write("physics.json", Physics(6, 4));
  Json values = Algorithm();
  values["MCAxisUpdate"] = true;
  values["Dbmps_max"] = 12;
  const auto algorithm = files.Write("algorithm.json", values);
  const EnhancedVMCUpdateParams vmc(physics.c_str(), algorithm.c_str());
  const EnhancedMCMeasureParams measure(physics.c_str(), algorithm.c_str());
  CHECK(vmc.axis_update_params.MCAxisUpdateDmax == 12);
  CHECK(measure.axis_update_params.MCAxisUpdateDmax == 12);
  CHECK(vmc.axis_update_params.dmax_from_dbmps_max && measure.axis_update_params.dmax_from_dbmps_max);
  CHECK(vmc.axis_update_params.MCAxisUpdateDmin == 1);
  CHECK(vmc.axis_update_params.MCAxisUpdateTruncErr == 0.0);
  CheckSliceParams(*Visit<Dense>(vmc.axis_update_params, vmc.mc_params).axis, 1, 12, 0.0);
}

void TestInvalidCombinationsThrow(const TempFiles &files) {
  const auto open = files.Write("physics.json", Physics(6, 4));
  // The std::invalid_argument messages of the VMC and the measurement parser.
  const auto messages = [&](const std::string &physics, const Json &values) {
    const auto algorithm = files.Write("algorithm.json", values);
    return std::pair<std::string, std::string>{
        InvalidArgumentMessage([&] { EnhancedVMCUpdateParams params(physics.c_str(), algorithm.c_str()); }),
        InvalidArgumentMessage([&] { EnhancedMCMeasureParams params(physics.c_str(), algorithm.c_str()); })};
  };
  // Both parsers reject an enabled axis update with message @p part.
  const auto rejects = [&](const std::string &physics, const Json &values, const std::string &part) {
    const auto [vmc, measure] = messages(physics, values);
    return Contains(vmc, part) && Contains(measure, part);
  };
  // Neither parser's message contains @p part (both must still throw).
  const auto omits = [&](const std::string &physics, const Json &values, const std::string &part) {
    const auto [vmc, measure] = messages(physics, values);
    return !vmc.empty() && !measure.empty() && !Contains(vmc, part) && !Contains(measure, part);
  };
  // Messages about MCAxisUpdateDmax say this only when the key is absent (default Dbmps_max).
  const std::string kKeyAbsent = "(the key is absent, so it takes its default Dbmps_max)";
  const auto accepts = [&](const std::string &physics, const Json &values) {
    const auto algorithm = files.Write("algorithm.json", values);
    EnhancedVMCUpdateParams vmc(physics.c_str(), algorithm.c_str());
    EnhancedMCMeasureParams measure(physics.c_str(), algorithm.c_str());
    return vmc.axis_update_params.MCAxisUpdate == measure.axis_update_params.MCAxisUpdate;
  };
  Json enabled = Algorithm();
  enabled["MCAxisUpdate"] = true;

  // Periodic boundaries.
  const auto periodic = files.Write("periodic.json", Physics(4, 4, "Periodic"));
  Json pbc = enabled;
  pbc.erase("Dbmps_max");
  pbc["MCAxisUpdateDmax"] = 8;
  CHECK(rejects(periodic, pbc, "requires BoundaryCondition=Open"));

  // Spin-inversion projection (the lattice satisfies the projection's own requirements).
  Json parity = enabled;
  parity["SpinInversionParity"] = 1;
  CHECK(rejects(open, parity, "not supported with SpinInversionParity != 0"));

  // MCAxisUpdateDmax = 0, explicit (with Dbmps_max = 10 present) or through an absent Dbmps_max.
  // Only the second message blames the default.
  Json zero = enabled;
  zero["MCAxisUpdateDmax"] = 0;
  CHECK(rejects(open, zero, "requires MCAxisUpdateDmax > 0, got MCAxisUpdateDmax = 0."));
  CHECK(omits(open, zero, kKeyAbsent));
  CHECK(omits(open, zero, "Dbmps_max"));
  Json no_bmps = enabled;
  no_bmps.erase("Dbmps_max");
  CHECK(rejects(open, no_bmps, "requires MCAxisUpdateDmax > 0, got MCAxisUpdateDmax = 0 " + kKeyAbsent +
                                   ". Set MCAxisUpdateDmax, or Dbmps_max."));

  // MCAxisUpdateDmin outside [1, MCAxisUpdateDmax], MCAxisUpdateTruncErr outside [0, 1).
  Json dmin = enabled;
  dmin["MCAxisUpdateDmin"] = 0;
  CHECK(rejects(open, dmin, "1 <= MCAxisUpdateDmin <= MCAxisUpdateDmax"));
  dmin["MCAxisUpdateDmin"] = 11;  // MCAxisUpdateDmax absent: Dbmps_max = 10
  CHECK(rejects(open, dmin, "got MCAxisUpdateDmin = 11 and MCAxisUpdateDmax = 10 " + kKeyAbsent + "."));
  dmin["MCAxisUpdateDmax"] = 10;  // the same value, set explicitly
  CHECK(rejects(open, dmin, "got MCAxisUpdateDmin = 11 and MCAxisUpdateDmax = 10."));
  dmin.erase("MCAxisUpdateDmax");
  dmin["MCAxisUpdateDmin"] = 10;
  CHECK(accepts(open, dmin));
  for (double error : {1.0, -1e-3}) {
    Json trunc = enabled;
    trunc["MCAxisUpdateTruncErr"] = error;
    CHECK(rejects(open, trunc, "MCAxisUpdateTruncErr must be in [0, 1)"));
  }

  // Below floor(max(Lx, Ly) / 2) + 1 with MCRestrictU1=true (the default): rows longer (Lx = 7)
  // and columns longer (Ly = 7) both give the bound 4.
  for (const auto &[lx, ly] : {std::pair<size_t, size_t>{7, 4}, std::pair<size_t, size_t>{4, 7}}) {
    const auto physics = files.Write("bound.json", Physics(lx, ly));
    Json bound = enabled;
    bound["MCAxisUpdateDmax"] = 3;
    CHECK(rejects(physics, bound, "MCAxisUpdateDmax = 3 is below floor(max(Lx, Ly) / 2) + 1 = 4"));
    CHECK(rejects(physics, bound, "section 5.4"));
    CHECK(rejects(physics, bound, "Set MCAxisUpdateDmax to at least 4."));
    CHECK(omits(physics, bound, "Dbmps_max"));  // Dbmps_max = 10 is present but not the cause
    bound["MCAxisUpdateDmax"] = 4;
    CHECK(accepts(physics, bound));

    // The same value through the default: Dbmps_max = 3 and no MCAxisUpdateDmax key.
    Json defaulted = enabled;
    defaulted["Dbmps_max"] = 3;
    CHECK(rejects(physics, defaulted,
                  "MCAxisUpdateDmax = 3 " + kKeyAbsent + " is below floor(max(Lx, Ly) / 2) + 1 = 4"));
    CHECK(rejects(physics, defaulted,
                  "Set MCAxisUpdateDmax to at least 4, or raise Dbmps_max, which supplies its default."));
    defaulted["MCAxisUpdateDmax"] = 4;  // an explicit key overrides the small Dbmps_max
    CHECK(accepts(physics, defaulted));
  }

  // The bound applies when the build's tensors carry S_z or MCRestrictU1=true, not otherwise.
  Json dense_free = enabled;
  dense_free["MCRestrictU1"] = false;
  dense_free["MCAxisUpdateDmax"] = 3;
  const auto dense_free_path = files.Write("dense_free.json", dense_free);
  const heisenberg_params::PhysicalParams physical(open.c_str());
  const heisenberg_params::MonteCarloNumericalParams mc(dense_free_path.c_str());
  const heisenberg_params::BMPSParams bmps(dense_free_path.c_str());
  CHECK(heisenberg_params::AxisUpdateParams::SectorBlockBound(physical) == 4);
  CHECK(InvalidArgumentMessage([&] {
          heisenberg_params::AxisUpdateParams axis(dense_free_path.c_str(), physical, mc, bmps, false, false);
        }).empty());
  CHECK(Contains(InvalidArgumentMessage([&] {
                   heisenberg_params::AxisUpdateParams axis(dense_free_path.c_str(), physical, mc, bmps, false,
                                                            /*u1_symmetric_build=*/true);
                 }),
                 "below floor(max(Lx, Ly) / 2) + 1 = 4"));
  // The parsers use this build's flag: a -DU1SYM build rejects the same file.
  if constexpr (heisenberg_params::kU1SymmetricBuild) {
    CHECK(rejects(open, dense_free, "below floor(max(Lx, Ly) / 2) + 1 = 4"));
  } else {
    CHECK(accepts(open, dense_free));
  }

  // Keys of the wrong JSON type are rejected even when the axis update is off.
  Json wrong = Algorithm();
  wrong["MCAxisUpdate"] = "true";
  CHECK(rejects(open, wrong, "'MCAxisUpdate' must be a boolean"));
  wrong = Algorithm();
  wrong["MCAxisUpdateDmax"] = -1;
  CHECK(rejects(open, wrong, "'MCAxisUpdateDmax' must be a non-negative integer"));
  wrong["MCAxisUpdateDmax"] = 2.5;
  CHECK(rejects(open, wrong, "'MCAxisUpdateDmax' must be a non-negative integer"));

  // Disabled, the values are not validated: no throw for PBC, Dmax = 0, Dmin > Dmax, bad error.
  Json disabled = pbc;
  disabled["MCAxisUpdate"] = false;
  disabled["MCAxisUpdateDmax"] = 0;
  disabled["MCAxisUpdateDmin"] = 5;
  disabled["MCAxisUpdateTruncErr"] = 2.0;
  disabled["SpinInversionParity"] = 0;
  CHECK(accepts(periodic, disabled));
}

}  // namespace axis_update_params_test

int main() {
  using namespace axis_update_params_test;
  try {
    const TempFiles files;
    TestDefaultsKeepTheLocalUpdater(files);
    TestExplicitValuesAndMapping(files);
    TestDmaxDefaultsToDbmpsMax(files);
    TestInvalidCombinationsThrow(files);
  } catch (const std::exception &error) {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  std::cout << "test_axis_update_params: all checks passed" << std::endl;
  return 0;
}
