// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/ond_dim_tn/boundary_mps/bmps.h"
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
    // IO defaults with overrides if provided
    wavefunction_base = this->ParseStrOr("WavefunctionBase", wavefunction_base);
    configuration_load_dir = this->ParseStrOr("ConfigurationLoadDir", configuration_load_dir);
    configuration_dump_dir = this->ParseStrOr("ConfigurationDumpDir", configuration_dump_dir);
    // Default config dirs to base+"final" (e.g., tpsfinal/)
    if (configuration_load_dir.empty()) configuration_load_dir = wavefunction_base + "final";
    if (configuration_dump_dir.empty()) configuration_dump_dir = wavefunction_base + "final";
  }

  heisenberg_params::PhysicalParams physical_params;        // Lx, Ly, J2, RemoveCorner, ModelType
  heisenberg_params::MonteCarloNumericalParams mc_params;   // MC_samples, WarmUp, MCLocalUpdateSweepsBetweenSample, MCRestrictU1
  heisenberg_params::BMPSParams bmps_params;                // Dbmps_min/max, MPSCompressScheme, TruncErr, ThreadNum

  // IO configuration (aligned with VMC)
  std::string wavefunction_base = "tps";        ///< basename for split-index TPS IO: base+"final", base+"lowest"
  std::string configuration_load_dir = "tpsfinal";      ///< directory to load MC configuration{rank}
  std::string configuration_dump_dir = "tpsfinal";      ///< directory to dump MC configuration{rank}

  /**
   * @brief Create PEPSParams (BMPS for OBC, TRG for PBC).
   */
  qlpeps::PEPSParams CreatePEPSParams() {
    if (physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic) {
      // Enforce TRG parameter presence for PBC runs.
      if (!(this->Has("TRGDmin") && this->Has("TRGDmax") && this->Has("TRGTruncErr"))) {
        throw std::invalid_argument(
            "PBC requested but TRG params are missing in algorithm JSON. "
            "Require: TRGDmin, TRGDmax, TRGTruncErr (optional: TRGInvRelativeEps).");
      }
      const size_t d_min = static_cast<size_t>(ParseInt("TRGDmin"));
      const size_t d_max = static_cast<size_t>(ParseInt("TRGDmax"));
      const double trunc_err = ParseDouble("TRGTruncErr");
      const double inv_eps = this->ParseDoubleOr("TRGInvRelativeEps", 1e-12);
      return qlpeps::PEPSParams(qlpeps::TRGTruncateParams<qlten::QLTEN_Double>(d_min, d_max, trunc_err, inv_eps));
    }
    if (!bmps_params.HasBMPSRequiredKeys()) {
      throw std::invalid_argument(
          "OBC requested but BMPS params are missing in algorithm JSON. "
          "Require: Dbmps_max (Dbmps_min optional; MPSCompressScheme optional, default=SVD).");
    }
    return qlpeps::PEPSParams(bmps_params.CreateTruncatePara());
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_MEASURE_PARAMS_PARSER_H
