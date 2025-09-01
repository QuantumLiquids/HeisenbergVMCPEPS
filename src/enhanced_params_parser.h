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
    
    // Parse optimizer configuration (backward-compatible defaults)
    optimizer_type = ParseStringDefault("OptimizerType", "StochasticReconfiguration");
    // Aliases: allow short names
    if (optimizer_type == "SR" || optimizer_type == "sr") {
      optimizer_type = "StochasticReconfiguration";
    }
    max_iterations = ParseIntDefault("MaxIterations", ParseIntDefault("UpdateNum", 10));
    learning_rate = ParseDoubleDefault("LearningRate", ParseDoubleDefault("StepLengthFirst", 0.01));
    
    // Parse convergence criteria (optional)
    energy_tolerance = ParseDoubleDefault("EnergyTolerance", 0.0);
    gradient_tolerance = ParseDoubleDefault("GradientTolerance", 0.0);
    plateau_patience = ParseIntDefault("PlateauPatience", max_iterations);
    
    // Parse algorithm-specific parameters
    if (optimizer_type == "SGD") {
      momentum = ParseDoubleDefault("Momentum", 0.0);
      nesterov = ParseBoolDefault("Nesterov", false);
      weight_decay = ParseDoubleDefault("WeightDecay", 0.0);
    } else if (optimizer_type == "Adam") {
      beta1 = ParseDoubleDefault("Beta1", 0.9);
      beta2 = ParseDoubleDefault("Beta2", 0.999);
      epsilon = ParseDoubleDefault("Epsilon", 1e-8);
      weight_decay = ParseDoubleDefault("WeightDecay", 0.0);
    } else if (optimizer_type == "AdaGrad") {
      epsilon = ParseDoubleDefault("Epsilon", 1e-8);
      initial_accumulator = ParseDoubleDefault("InitialAccumulator", 0.0);
    } else if (optimizer_type == "StochasticReconfiguration") {
      // SR: require explicit CG params (no defaults) except Diagonal shift (defaults to 0.0)
      cg_max_iter = ParseInt("CGMaxIter");
      cg_tol = ParseDouble("CGTol");
      cg_residue_restart = ParseInt("CGResidueRestart");
      cg_diag_shift = ParseDoubleDefault("CGDiagShift", 0.0);
      normalize_update = ParseBoolDefault("NormalizeUpdate", false);
    }
    
    // Parse learning rate scheduler (optional)
    lr_scheduler_type = ParseStringDefault("LRScheduler", "");
    if (lr_scheduler_type == "ExponentialDecay") {
      decay_rate = ParseDouble("DecayRate");
      decay_steps = ParseInt("DecaySteps");
    } else if (lr_scheduler_type == "CosineAnnealing") {
      min_learning_rate = ParseDoubleDefault("MinLearningRate", 0.0);
    } else if (lr_scheduler_type == "Plateau") {
      plateau_factor = ParseDoubleDefault("PlateauFactor", 0.5);
      plateau_threshold = ParseDoubleDefault("PlateauThreshold", 1e-4);
    }
    
    // Parse gradient clipping (optional)
    clip_norm = ParseDoubleOptional("ClipNorm");
    clip_value = ParseDoubleOptional("ClipValue");
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
  
  /**
   * @brief Create qlpeps::VMCPEPSOptimizerParams with modern optimizer support
   */
  qlpeps::VMCPEPSOptimizerParams CreateVMCOptimizerParams(int rank = 0) {
    wavefunction_base = ParseStringDefault("WavefunctionBase", wavefunction_base);
    configuration_load_dir = ParseStringDefault("ConfigurationLoadDir", configuration_load_dir);
    configuration_dump_dir = ParseStringDefault("ConfigurationDumpDir", configuration_dump_dir);
    // Default config dirs to base+"final" (e.g., tpsfinal/)
    if (configuration_load_dir.empty()) configuration_load_dir = wavefunction_base + "final";
    if (configuration_dump_dir.empty()) configuration_dump_dir = wavefunction_base + "final";

    // Create Monte Carlo parameters
    size_t N = physical_params.Lx * physical_params.Ly;
    qlpeps::Configuration config(physical_params.Ly, physical_params.Lx);
    size_t spin_up_sites = N / 2;
    for (size_t row = 0; row < physical_params.Ly; row++) {
      for (size_t col = 0; col < physical_params.Lx; col++) {
        size_t site_idx = row * physical_params.Lx + col;
        config({row, col}) = (site_idx < spin_up_sites) ? 1 : 0;
      }
    }
    bool warmed_up = false;
    // Try load configuration from directory if provided: configuration{rank}
    if (!configuration_load_dir.empty()) {
      std::string load_dir = configuration_load_dir;
      if (!load_dir.empty() && load_dir.back() != '/') load_dir.push_back('/');
      std::string cfg_file = load_dir + "configuration" + std::to_string(rank);
      if (config.Load(cfg_file, 0)) {
        warmed_up = true;
      }
    }
    
    qlpeps::MonteCarloParams mc_params_obj(
      mc_params.MC_samples,
      mc_params.WarmUp,
      mc_params.MCLocalUpdateSweepsBetweenSample,
      config,
      warmed_up,
      configuration_dump_dir
    );
    
    qlpeps::PEPSParams peps_params_obj(
      bmps_params.CreateTruncatePara()
    );
    
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
    if (optimizer_type == "SGD") {
      qlpeps::SGDParams sgd_params(momentum, nesterov, weight_decay);
      return qlpeps::OptimizerParams(base_params, sgd_params);
    } else if (optimizer_type == "Adam") {
      qlpeps::AdamParams adam_params(beta1, beta2, epsilon, weight_decay);
      return qlpeps::OptimizerParams(base_params, adam_params);
    } else if (optimizer_type == "AdaGrad") {
      qlpeps::AdaGradParams adagrad_params(epsilon, initial_accumulator);
      return qlpeps::OptimizerParams(base_params, adagrad_params);
    } else if (optimizer_type == "StochasticReconfiguration") {
      qlpeps::ConjugateGradientParams cg_params(cg_max_iter, cg_tol, cg_residue_restart, cg_diag_shift);
      qlpeps::StochasticReconfigurationParams sr_params(cg_params, normalize_update);
      return qlpeps::OptimizerParams(base_params, sr_params);
    } else {
      throw std::invalid_argument("Unknown optimizer type: " + optimizer_type);
    }
  }
  
  // Helper methods for optional parameter parsing
  std::string ParseStringDefault(const std::string& key, const std::string& default_val) {
    try {
      return ParseStr(key);
    } catch (...) {
      return default_val;
    }
  }
  
  double ParseDoubleDefault(const std::string& key, double default_val) {
    try {
      return ParseDouble(key);
    } catch (...) {
      return default_val;
    }
  }
  
  size_t ParseIntDefault(const std::string& key, size_t default_val) {
    try {
      return ParseInt(key);
    } catch (...) {
      return default_val;
    }
  }
  
  bool ParseBoolDefault(const std::string& key, bool default_val) {
    try {
      return ParseBool(key);
    } catch (...) {
      return default_val;
    }
  }
  
  std::optional<double> ParseDoubleOptional(const std::string& key) {
    try {
      return ParseDouble(key);
    } catch (...) {
      return std::nullopt;
    }
  }
};

#endif // HEISENBERGVMCPEPS_ENHANCED_PARAMS_PARSER_H
