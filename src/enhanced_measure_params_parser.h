// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/ond_dim_tn/boundary_mps/bmps.h"
#include "common_params.h"

/**
 * @brief Enhanced Measure parameters that align with new two-file system.
 */
struct EnhancedMCMeasureParams {
  EnhancedMCMeasureParams(const char *physics_file, const char *algorithm_file)
      : physical_params(physics_file),
        mc_params(algorithm_file),
        bmps_params(algorithm_file) {}

  heisenberg_params::PhysicalParams physical_params;        // Lx, Ly, J2, RemoveCorner, ModelType, MCRestrictU1
  heisenberg_params::MonteCarloNumericalParams mc_params;   // MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample
  heisenberg_params::BMPSParams bmps_params;                // Dbmps_min/max, MPSCompressScheme, TruncErr, ThreadNum

  /**
   * @brief Create truncate para for BMPS using current numeric settings.
   */
  qlpeps::BMPSTruncatePara CreateBMPSPara() const {
    return bmps_params.CreateTruncatePara();
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H


