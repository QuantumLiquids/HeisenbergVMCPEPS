// SPDX-License-Identifier: LGPL-3.0-only
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "../src/enhanced_params_parser.h"
#include "../src/enhanced_measure_params_parser.h"
#include "../src/qldouble.h"

#ifdef USE_COMPLEX
static_assert(std::is_same_v<TenElemT, qlten::QLTEN_Complex>);
#else
static_assert(std::is_same_v<TenElemT, qlten::QLTEN_Double>);
#endif

namespace {
using qlpeps::config::Json;

void Check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

/** @brief Exercise both application parsers with the same physics and algorithm files. */
class Inputs {
 public:
  Inputs() : directory_(std::filesystem::temp_directory_path() /
      ("projection-params-" + std::to_string(
          std::chrono::steady_clock::now().time_since_epoch().count()))) {
    std::filesystem::create_directories(directory_);
  }
  ~Inputs() { std::filesystem::remove_all(directory_); }

  void Write(const Json &physics, const Json &algorithm) const {
    std::ofstream(Physics()) << Json{{"CaseParams", physics}};
    std::ofstream(Algorithm()) << Json{{"CaseParams", algorithm}};
  }
  std::string Physics() const { return (directory_ / "physics.json").string(); }
  std::string Algorithm() const { return (directory_ / "algorithm.json").string(); }

  void Accept(const qlpeps::PointGroupProjectionParams &expected) const {
    EnhancedVMCUpdateParams vmc(Physics().c_str(), Algorithm().c_str());
    const EnhancedMCMeasureParams measure(Physics().c_str(), Algorithm().c_str());
    Check(vmc.point_group_projection == expected, "Unexpected VMC projection");
    Check(measure.point_group_projection == expected, "VMC/measure projector mismatch");
    Check(vmc.spin_inversion_parity == expected.spin_inversion_parity, "Legacy parity mismatch");
    if (expected.group != "None") {
      const auto params = vmc.CreateVMCOptimizerParams();
      Check(params.mc_params.point_group_projection == expected, "Projection not passed to optimizer");
      Check(!params.mc_params.assume_initial_config_thermalized, "Projection skipped warm-up");
    }
  }
  void Reject() const {
    bool vmc_rejected = false, measure_rejected = false;
    try { EnhancedVMCUpdateParams params(Physics().c_str(), Algorithm().c_str()); }
    catch (const std::invalid_argument &) { vmc_rejected = true; }
    try { EnhancedMCMeasureParams params(Physics().c_str(), Algorithm().c_str()); }
    catch (const std::invalid_argument &) { measure_rejected = true; }
    Check(vmc_rejected && measure_rejected, "Invalid projection accepted by a driver");
  }

 private:
  std::filesystem::path directory_;
};
}  // namespace

int main() {
  Inputs files;
  const Json physics{{"Lx", 4}, {"Ly", 4}, {"J2", 0.0},
                     {"ModelType", "SquareHeisenberg"}, {"BoundaryCondition", "Open"}};
  const Json plain{{"MC_total_samples", 8}, {"WarmUp", 10},
                   {"MCLocalUpdateSweepsBetweenSample", 1}, {"Dbmps_max", 8},
                   {"OptimizerType", "SGD"}, {"ConfigurationLoadDir", "/nonexistent-projection-config"}};
  files.Write(physics, plain);
  files.Accept({});
  for (const std::string group : {"C4", "D4"}) {
    const std::vector<std::string> irreps = group == "C4"
        ? std::vector<std::string>{"0", "1", "2", "3"}
        : std::vector<std::string>{"A1", "A2", "B1", "B2", "E"};
    for (const auto &irrep : irreps) {
      for (int parity : {0, -1, 1}) {
        Json algorithm = plain;
        algorithm["PointGroup"] = group;
        algorithm["PointGroupIrrep"] = irrep;
        algorithm["SpinInversionParity"] = parity;
        files.Write(physics, algorithm);
#ifndef USE_COMPLEX
        if (group == "C4" && (irrep == "1" || irrep == "3")) {
          files.Reject();
          continue;
        }
#endif
        files.Accept({group, irrep, parity});
      }
    }
  }
  Json projected = plain;
  projected["PointGroup"] = "D4";
  files.Write(physics, projected);
  files.Accept({"D4", "A1", 0});
  Json c4 = plain;
  c4["PointGroup"] = "C4";
  files.Write(physics, c4);
  files.Accept({"C4", "0", 0});
  Json spin = plain;
  spin["SpinInversionParity"] = -1;
  files.Write(physics, spin);
  files.Accept({"None", "A1", -1});

  for (const auto &change : std::vector<Json>{
      {{"PointGroup", "D4h"}}, {{"PointGroup", 4}}, {{"PointGroupIrrep", "1"}},
      {{"PointGroupIrrep", 1}}, {{"SpinInversionParity", 2}},
      {{"SpinInversionParity", 0.5}}, {{"MCAxisUpdate", true}}, {{"MCRestrictU1", false}}}) {
    Json algorithm = projected;
    algorithm.update(change);
    files.Write(physics, algorithm);
    files.Reject();
  }
  for (const auto &change : std::vector<Json>{
      {{"Lx", 6}}, {{"Lx", 1}, {"Ly", 1}}, {{"BoundaryCondition", "Periodic"}},
      {{"RemoveCorner", true}}, {{"ModelType", "TriangleHeisenberg"}}}) {
    Json invalid_physics = physics;
    invalid_physics.update(change);
    files.Write(invalid_physics, projected);
    files.Reject();
  }
  Json odd = physics;
  odd["Lx"] = 3;
  odd["Ly"] = 3;
  files.Write(odd, projected);
  files.Accept({"D4", "A1", 0});
  projected["SpinInversionParity"] = 1;
  files.Write(odd, projected);
  files.Reject();
  Json unrestrained = projected;
  unrestrained["MCRestrictU1"] = false;
  files.Write(physics, unrestrained);
  files.Reject();
  for (const std::string model : {"SquareHeisenberg", "SquareXY"}) {
    Json j2_physics = physics;
    j2_physics["ModelType"] = model;
    j2_physics["J2"] = 0.5;
    files.Write(j2_physics, projected);
    files.Accept({"D4", "A1", 1});
    files.Write(j2_physics, spin);
    files.Reject();
  }
  // Disabled projection preserves the rectangular unprojected model path.
  Json rectangle = physics;
  rectangle["Lx"] = 6;
  files.Write(rectangle, plain);
  files.Accept({});
}
