//
// Enhanced parameter parser with modern optimizer support
//

#ifndef HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer_params.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "qlpeps/optimizer/optimizer_params.h"
#include "qlpeps/api/config/optimizer_params_parser.h"
#include "common_params.h"
#include "projection_params.h"
#include <algorithm>
#include <cctype>
#include <optional>
#include <fstream>

namespace heisenberg_params {
/// Keep application defaults and temporary MinSR aliases outside the public parser.
inline qlpeps::OptimizerParams ReadOptimizerParams(const char *file) {
  auto params = qlpeps::config::ReadCaseParamsFile(file);
  if (!params.contains("OptimizerType")) params["OptimizerType"] = "StochasticReconfiguration";
  if (!params.contains("MaxIterations")) params["MaxIterations"] = 10;
  if (!params.contains("LearningRate")) params["LearningRate"] = 0.01;
  qlpeps::config::PromoteParameterAlias(params, "MinSRRelativePInv", "MinSRRPinv");
  qlpeps::config::PromoteParameterAlias(params, "MinSRAbsolutePInv", "MinSRAPinv");
  return qlpeps::config::ParseOptimizerParams(params);
}
}  // namespace heisenberg_params

/**
 * @brief Model inputs composed with the shared PEPS optimizer configuration.
 */
struct EnhancedVMCUpdateParams : public qlmps::CaseParamsParserBasic {
  EnhancedVMCUpdateParams(const char *physics_file, const char *algorithm_file) : 
      CaseParamsParserBasic(algorithm_file),
      physical_params(physics_file),
      mc_params(algorithm_file),
      bmps_params(algorithm_file),
      optimizer_params(heisenberg_params::ReadOptimizerParams(algorithm_file)) {
    point_group_projection = heisenberg_params::ReadProjectionParams(
        algorithm_file, physical_params, mc_params);
    spin_inversion_parity = point_group_projection.spin_inversion_parity;
    io_params.Parse(*this);
    axis_update_params = heisenberg_params::AxisUpdateParams(
        algorithm_file, physical_params, mc_params, bmps_params, spin_inversion_parity != 0);
    heisenberg_params::ValidatePointGroupUpdater(point_group_projection, axis_update_params);
  }

  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  
  /// Zero preserves plain PEPS; +/-1 selects psi(x) +/- psi(Fx).
  int spin_inversion_parity = 0;

  /// Spatial group/irrep and independent spin-inversion parity.
  qlpeps::PointGroupProjectionParams point_group_projection;

  /// Opt-in axis update of OBC sampling (`MCAxisUpdate*` keys); disabled by default.
  heisenberg_params::AxisUpdateParams axis_update_params;

  qlpeps::OptimizerParams optimizer_params;
  heisenberg_params::IOParams io_params;

  /**
   * @brief Create qlpeps::VMCPEPSOptimizerParams with modern optimizer support
   */
  qlpeps::VMCPEPSOptimizerParams CreateVMCOptimizerParams(int rank = 0) {
    // Create Monte Carlo parameters
    auto [config, warmed_up] = heisenberg_params::InitOrLoadConfigWithStrategy(
        physical_params,
        mc_params,
        io_params.configuration_load_dir,
        rank);

    auto mc_params_obj = mc_params.CreateMonteCarloParams(
        config, warmed_up, io_params.configuration_dump_dir);
    if (point_group_projection.group != "None") {
      mc_params_obj.point_group_projection = point_group_projection;
      mc_params_obj.assume_initial_config_thermalized = false;
    }

    const qlpeps::ContractorParams contractor_params_obj = heisenberg_params::CreateContractorParams(
        physical_params.BoundaryCondition, bmps_params, bmps_params.algorithm_values);

    return qlpeps::VMCPEPSOptimizerParams(optimizer_params, mc_params_obj, contractor_params_obj);
  }

};

#endif // HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H
