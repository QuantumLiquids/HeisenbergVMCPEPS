// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/one_dim_tn/boundary_mps/bmps_truncate_params.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "common_params.h"
#include "projection_params.h"

/**
 * @brief Enhanced Measure parameters that align with new two-file system.
 */
struct EnhancedMCMeasureParams : public qlmps::CaseParamsParserBasic {
  EnhancedMCMeasureParams(const char *physics_file, const char *algorithm_file)
      : qlmps::CaseParamsParserBasic(algorithm_file),
        physical_params(physics_file),
        mc_params(algorithm_file),
        bmps_params(algorithm_file) {
    io_params.Parse(*this);
    point_group_projection = heisenberg_params::ReadProjectionParams(
        algorithm_file, physical_params, mc_params);
    axis_update_params = heisenberg_params::AxisUpdateParams(
        algorithm_file, physical_params, mc_params, bmps_params,
        point_group_projection.spin_inversion_parity != 0);
    heisenberg_params::ValidatePointGroupUpdater(point_group_projection, axis_update_params);
  }

  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  heisenberg_params::IOParams io_params;
  /// Same projector used for VMC and all measured local estimators.
  qlpeps::PointGroupProjectionParams point_group_projection;
  /// Opt-in axis update of OBC sampling (`MCAxisUpdate*` keys); disabled by default.
  heisenberg_params::AxisUpdateParams axis_update_params;

  /**
   * @brief Create ContractorParams (BMPS for OBC, TRG or HOTRG for PBC).
   */
  qlpeps::ContractorParams CreateContractorParams() {
    return heisenberg_params::CreateContractorParams(
        physical_params.BoundaryCondition, bmps_params, bmps_params.algorithm_values);
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
