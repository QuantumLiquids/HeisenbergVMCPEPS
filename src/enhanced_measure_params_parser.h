// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/one_dim_tn/boundary_mps/bmps.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "common_params.h"

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
  }

  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  heisenberg_params::IOParams io_params;

  /**
   * @brief Create PEPSParams (BMPS for OBC, TRG for PBC).
   */
  qlpeps::PEPSParams CreatePEPSParams() {
    return heisenberg_params::CreatePEPSParams(
        physical_params.BoundaryCondition, bmps_params, *this);
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
