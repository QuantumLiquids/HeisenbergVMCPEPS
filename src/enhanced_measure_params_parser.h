// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/ond_dim_tn/boundary_mps/bmps.h"
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
    // IO defaults with overrides if provided
    wavefunction_base = this->ParseStrOr("WavefunctionBase", wavefunction_base);
    configuration_load_dir = this->ParseStrOr("ConfigurationLoadDir", configuration_load_dir);
    configuration_dump_dir = this->ParseStrOr("ConfigurationDumpDir", configuration_dump_dir);
    // Default config dirs to base+"final" (e.g., tpsfinal/)
    if (configuration_load_dir.empty()) configuration_load_dir = wavefunction_base + "final";
    if (configuration_dump_dir.empty()) configuration_dump_dir = wavefunction_base + "final";
  }

  heisenberg_params::PhysicalParams physical_params;        // Lx, Ly, J2, RemoveCorner, ModelType, MCRestrictU1
  heisenberg_params::MonteCarloNumericalParams mc_params;   // MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample
  heisenberg_params::BMPSParams bmps_params;                // Dbmps_min/max, MPSCompressScheme, TruncErr, ThreadNum

  // IO configuration (aligned with VMC)
  std::string wavefunction_base = "tps";        ///< basename for split-index TPS IO: base+"final", base+"lowest"
  std::string configuration_load_dir = "tpsfinal";      ///< directory to load MC configuration{rank}
  std::string configuration_dump_dir = "tpsfinal";      ///< directory to dump MC configuration{rank}

  /**
   * @brief Create truncate para for BMPS using current numeric settings.
   */
  qlpeps::BMPSTruncatePara CreateBMPSPara() const {
    return bmps_params.CreateTruncatePara();
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H


