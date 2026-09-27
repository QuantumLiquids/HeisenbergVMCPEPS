// SPDX-License-Identifier: LGPL-3.0-only
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include "../src/enhanced_params_parser.h"

void Require(bool condition) {
  if (!condition) throw std::runtime_error("Optimizer adapter regression");
}

int main() {
  using qlpeps::config::Json;
  const auto file = std::filesystem::temp_directory_path() /
      ("optimizer-adapter-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
  struct Cleanup { std::filesystem::path path; ~Cleanup() { std::filesystem::remove(path); } } cleanup{file};
  Json input{{"OptimizerType", "MinSR"}, {"MaxIterations", 17}, {"LearningRate", 0.04},
             {"MinSRRelativePInv", 5e-9}, {"MinSRAbsolutePInv", 3e-10}};
  auto read = [&] {
    std::ofstream out(file);
    out << Json{{"CaseParams", input}};
    out.close();
    return heisenberg_params::ReadOptimizerParams(file.c_str());
  };
  auto rejects = [&] {
    try { (void) read(); } catch (const std::invalid_argument &) { return true; }
    return false;
  };
  const auto canonical = read();
  Require(canonical.GetAlgorithmParams<qlpeps::MinSRParams>().r_pinv == 5e-9);
  Require(canonical.GetAlgorithmParams<qlpeps::MinSRParams>().a_pinv == 3e-10);
  Require(canonical.base_params.max_iterations == 17);
  Require(canonical.base_params.learning_rate == 0.04);
  input.erase("MinSRRelativePInv");
  input.erase("MinSRAbsolutePInv");
  input["MinSRRPinv"] = 5e-9;
  input["MinSRAPinv"] = 3e-10;
  const auto legacy = read();
  Require(legacy.GetAlgorithmParams<qlpeps::MinSRParams>().r_pinv == 5e-9);
  Require(legacy.GetAlgorithmParams<qlpeps::MinSRParams>().a_pinv == 3e-10);
  input["MinSRRelativePInv"] = 5e-9;
  Require(rejects());
  input = Json{{"OptimizerType", "SGD"}};
  const auto defaults = read();
  Require(defaults.base_params.max_iterations == 10);
  Require(defaults.base_params.learning_rate == 0.01);
}
