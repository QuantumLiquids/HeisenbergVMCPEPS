//
// Enhanced parameter parser with modern optimizer support
//

#ifndef HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer_params.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "qlpeps/optimizer/optimizer_params.h"
#include "qlpeps/optimizer/lr_schedulers.h"
#include "common_params.h"
#include <optional>
#include <fstream>

/**
 * @brief Enhanced VMC parameters with modern optimizer support
 */
struct EnhancedVMCUpdateParams : public qlmps::CaseParamsParserBasic {
  EnhancedVMCUpdateParams(const char *physics_file, const char *algorithm_file) : 
      CaseParamsParserBasic(algorithm_file),
      physical_params(physics_file),
      mc_params(algorithm_file),
      bmps_params(algorithm_file) {
    
    // Parse optimizer configuration with default values
    optimizer_type = this->ParseStrOr("OptimizerType", "StochasticReconfiguration");
    // Aliases: allow short names
    if (optimizer_type == "SR" || optimizer_type == "sr") {
      optimizer_type = "StochasticReconfiguration";
    }
    max_iterations = static_cast<size_t>(this->ParseIntOr("MaxIterations", 10));
    learning_rate = this->ParseDoubleOr("LearningRate", 0.01);
    
    // Parse convergence criteria (optional)
    // Convergence criteria (optional)
    energy_tolerance = this->ParseDoubleOr("EnergyTolerance", 0.0);
    gradient_tolerance = this->ParseDoubleOr("GradientTolerance", 0.0);
    plateau_patience = static_cast<size_t>(this->ParseIntOr("PlateauPatience", static_cast<int>(max_iterations)));
    
    // Parse algorithm-specific parameters
    if (optimizer_type == "SGD") {
      momentum = this->ParseDoubleOr("Momentum", 0.0);
      nesterov = this->ParseBoolOr("Nesterov", false);
      weight_decay = this->ParseDoubleOr("WeightDecay", 0.0);
    } else if (optimizer_type == "Adam") {
      beta1 = this->ParseDoubleOr("Beta1", 0.9);
      beta2 = this->ParseDoubleOr("Beta2", 0.999);
      epsilon = this->ParseDoubleOr("Epsilon", 1e-8);
      weight_decay = this->ParseDoubleOr("WeightDecay", 0.0);
    } else if (optimizer_type == "AdaGrad") {
      epsilon = this->ParseDoubleOr("Epsilon", 1e-8);
      initial_accumulator = this->ParseDoubleOr("InitialAccumulator", 0.0);
    } else if (optimizer_type == "StochasticReconfiguration") {
      // SR: require explicit CG params (no defaults) except Diagonal shift (defaults to 0.0)
      cg_max_iter = ParseInt("CGMaxIter");
      cg_tol = ParseDouble("CGTol");
      cg_residue_restart = ParseInt("CGResidueRestart");
      cg_diag_shift = this->ParseDoubleOr("CGDiagShift", 0.0);
      normalize_update = this->ParseBoolOr("NormalizeUpdate", false);
    }
    
    // Parse learning rate scheduler (optional)
    // Learning rate scheduler (optional)
    lr_scheduler_type = this->ParseStrOr("LRScheduler", "");
    if (lr_scheduler_type == "ExponentialDecay") {
      decay_rate = ParseDouble("DecayRate");
      decay_steps = ParseInt("DecaySteps");
    } else if (lr_scheduler_type == "CosineAnnealing") {
      min_learning_rate = this->ParseDoubleOr("MinLearningRate", 0.0);
    } else if (lr_scheduler_type == "Plateau") {
      plateau_factor = this->ParseDoubleOr("PlateauFactor", 0.5);
      plateau_threshold = this->ParseDoubleOr("PlateauThreshold", 1e-4);
    }
    
    // Parse gradient clipping (optional)
    // Gradient clipping (optional)
    double tmp_val = 0.0;
    if (this->TryParseDouble("ClipNorm", tmp_val)) clip_norm = tmp_val; else clip_norm.reset();
    if (this->TryParseDouble("ClipValue", tmp_val)) clip_value = tmp_val; else clip_value.reset();

    // Spike recovery (optional; default: disabled to preserve prior behavior)
    spike_auto_recover = this->ParseBoolOr("SpikeAutoRecover", false);
    spike_max_retries = static_cast<size_t>(this->ParseIntOr("SpikeMaxRetries", 2));
    spike_factor_err = this->ParseDoubleOr("SpikeFactorErr", 100.0);
    spike_factor_grad = this->ParseDoubleOr("SpikeFactorGrad", 1e10);
    spike_factor_ngrad = this->ParseDoubleOr("SpikeFactorNGrad", 10.0);
    spike_sr_min_iters_suspicious = static_cast<size_t>(this->ParseIntOr("SpikeSRMinItersSuspicious", 1));
    spike_enable_rollback = this->ParseBoolOr("SpikeEnableRollback", false);
    spike_ema_window = static_cast<size_t>(this->ParseIntOr("SpikeEMAWindow", 50));
    spike_sigma_k = this->ParseDoubleOr("SpikeSigmaK", 10.0);
    spike_log_csv = this->ParseStrOr("SpikeLogCSV", "");
  }

  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  
  // Optimizer configuration
  std::string optimizer_type;
  size_t max_iterations;
  double learning_rate;
  double energy_tolerance;
  double gradient_tolerance;
  size_t plateau_patience;
  
  // IO configuration
  std::string wavefunction_base = "tps";        ///< basename for split-index TPS IO: base+"final", base+"lowest"
  std::string configuration_load_dir = "tpsfinal";      ///< directory to load MC configuration{rank}; 
  std::string configuration_dump_dir = "tpsfinal";      ///< directory to dump MC configuration{rank}; 
  
  // SGD parameters
  double momentum = 0.0;
  bool nesterov = false;
  double weight_decay = 0.0;
  
  // Adam parameters
  double beta1 = 0.9;
  double beta2 = 0.999;
  double epsilon = 1e-8;
  
  // AdaGrad parameters
  double initial_accumulator = 0.0;
  
  // Stochastic Reconfiguration parameters
  size_t cg_max_iter = 100;
  double cg_tol = 1e-8;
  size_t cg_residue_restart = 20;
  double cg_diag_shift = 0.01;
  bool normalize_update = false;
  
  // Learning rate scheduler
  std::string lr_scheduler_type;
  double decay_rate = 0.95;
  size_t decay_steps = 10;
  double min_learning_rate = 0.0;
  double plateau_factor = 0.5;
  double plateau_threshold = 1e-4;
  
  // Gradient clipping
  std::optional<double> clip_norm;
  std::optional<double> clip_value;

  // Spike recovery
  bool spike_auto_recover = false;
  size_t spike_max_retries = 2;
  double spike_factor_err = 100.0;
  double spike_factor_grad = 1e10;
  double spike_factor_ngrad = 10.0;
  size_t spike_sr_min_iters_suspicious = 1;
  bool spike_enable_rollback = false;
  size_t spike_ema_window = 50;
  double spike_sigma_k = 10.0;
  std::string spike_log_csv;

  /**
   * @brief Create qlpeps::VMCPEPSOptimizerParams with modern optimizer support
   */
  qlpeps::VMCPEPSOptimizerParams CreateVMCOptimizerParams(int rank = 0) {
    // IO defaults with overrides if provided
    wavefunction_base = this->ParseStrOr("WavefunctionBase", wavefunction_base);
    configuration_load_dir = this->ParseStrOr("ConfigurationLoadDir", configuration_load_dir);
    configuration_dump_dir = this->ParseStrOr("ConfigurationDumpDir", configuration_dump_dir);
    // Default config dirs to base+"final" (e.g., tpsfinal/)
    if (configuration_load_dir.empty()) configuration_load_dir = wavefunction_base + "final";
    if (configuration_dump_dir.empty()) configuration_dump_dir = wavefunction_base + "final";

    // Create Monte Carlo parameters
    size_t N = physical_params.Lx * physical_params.Ly;
    size_t spin_up_sites = N / 2;
    qlpeps::OccupancyNum occupancy = {spin_up_sites, N - spin_up_sites};
    qlpeps::Configuration config(physical_params.Ly, physical_params.Lx, occupancy);
    bool warmed_up = false;
    // Try load configuration from directory if provided: configuration{rank}
    if (!configuration_load_dir.empty()) {
      warmed_up = config.Load(configuration_load_dir, static_cast<size_t>(rank));
    }
    
    qlpeps::MonteCarloParams mc_params_obj(
      mc_params.MC_total_samples,
      mc_params.WarmUp,
      mc_params.MCLocalUpdateSweepsBetweenSample,
      config,
      warmed_up,
      configuration_dump_dir
    );
    
    // Backend selection:
    // - OBC uses BMPS truncation parameters.
    // - PBC uses TRG truncation parameters.
    const qlpeps::PEPSParams peps_params_obj = [&]() -> qlpeps::PEPSParams {
      if (physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic) {
        // Enforce TRG parameter presence for PBC runs. Do not silently fall back to BMPS knobs.
        if (!(this->Has("TRGDmin") && this->Has("TRGDmax") && this->Has("TRGTruncErr"))) {
          throw std::invalid_argument(
              "PBC requested but TRG params are missing in algorithm JSON. "
              "Require: TRGDmin, TRGDmax, TRGTruncErr (optional: TRGInvRelativeEps).");
        }
        const size_t d_min = static_cast<size_t>(ParseInt("TRGDmin"));
        const size_t d_max = static_cast<size_t>(ParseInt("TRGDmax"));
        const double trunc_err = ParseDouble("TRGTruncErr");
        const double inv_eps = this->ParseDoubleOr("TRGInvRelativeEps", 1e-12);
        return qlpeps::PEPSParams(
            qlpeps::TRGTruncateParams<qlten::QLTEN_Double>(d_min, d_max, trunc_err, inv_eps));
      }

      if (!bmps_params.HasBMPSRequiredKeys()) {
        throw std::invalid_argument(
            "OBC requested but BMPS params are missing in algorithm JSON. "
            "Require: Dbmps_max (Dbmps_min optional; MPSCompressScheme optional, default=SVD).");
      }
      return qlpeps::PEPSParams(bmps_params.CreateTruncatePara());
    }();
    
    // Create optimizer parameters based on type
    qlpeps::OptimizerParams opt_params = CreateOptimizerParams();
    
    return qlpeps::VMCPEPSOptimizerParams(opt_params, mc_params_obj, peps_params_obj);
  }

private:
  /**
   * @brief Create type-specific optimizer parameters
   */
  qlpeps::OptimizerParams CreateOptimizerParams() {
    // Create learning rate scheduler if specified
    std::unique_ptr<qlpeps::LearningRateScheduler> scheduler = nullptr;
    if (lr_scheduler_type == "ExponentialDecay") {
      scheduler = std::make_unique<qlpeps::ExponentialDecayLR>(learning_rate, decay_rate, decay_steps);
    } else if (lr_scheduler_type == "CosineAnnealing") {
      scheduler = std::make_unique<qlpeps::CosineAnnealingLR>(learning_rate, max_iterations, min_learning_rate);
    } else if (lr_scheduler_type == "Plateau") {
      scheduler = std::make_unique<qlpeps::PlateauLR>(learning_rate, plateau_factor, plateau_patience, plateau_threshold);
    }
    
    // Create base parameters
    qlpeps::OptimizerParams::BaseParams base_params(
      max_iterations, energy_tolerance, gradient_tolerance, plateau_patience, 
      learning_rate, std::move(scheduler)
    );
    
    // Set gradient clipping if specified
    if (clip_norm) base_params.clip_norm = *clip_norm;
    if (clip_value) base_params.clip_value = *clip_value;
    
    // Create algorithm-specific parameters
    qlpeps::CheckpointParams ckpt_params{};
    qlpeps::SpikeRecoveryParams spike_params = CreateSpikeRecoveryParams();

    if (optimizer_type == "SGD") {
      qlpeps::SGDParams sgd_params(momentum, nesterov, weight_decay);
      return qlpeps::OptimizerParams(base_params, sgd_params, ckpt_params, spike_params);
    } else if (optimizer_type == "Adam") {
      qlpeps::AdamParams adam_params(beta1, beta2, epsilon, weight_decay);
      return qlpeps::OptimizerParams(base_params, adam_params, ckpt_params, spike_params);
    } else if (optimizer_type == "AdaGrad") {
      qlpeps::AdaGradParams adagrad_params(epsilon, initial_accumulator);
      return qlpeps::OptimizerParams(base_params, adagrad_params, ckpt_params, spike_params);
    } else if (optimizer_type == "StochasticReconfiguration") {
      qlpeps::ConjugateGradientParams cg_params(cg_max_iter, cg_tol, cg_residue_restart, cg_diag_shift);
      qlpeps::StochasticReconfigurationParams sr_params(cg_params, normalize_update);
      return qlpeps::OptimizerParams(base_params, sr_params, ckpt_params, spike_params);
    } else {
      throw std::invalid_argument("Unknown optimizer type: " + optimizer_type);
    }
  }

  qlpeps::SpikeRecoveryParams CreateSpikeRecoveryParams() {
    qlpeps::SpikeRecoveryParams p;
    p.enable_auto_recover = spike_auto_recover;
    p.redo_mc_max_retries = spike_max_retries;
    p.factor_err = spike_factor_err;
    p.factor_grad = spike_factor_grad;
    p.factor_ngrad = spike_factor_ngrad;
    p.sr_min_iters_suspicious = spike_sr_min_iters_suspicious;
    p.enable_rollback = spike_enable_rollback;
    p.ema_window = spike_ema_window;
    p.sigma_k = spike_sigma_k;
    p.log_trigger_csv_path = spike_log_csv;
    return p;
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H
