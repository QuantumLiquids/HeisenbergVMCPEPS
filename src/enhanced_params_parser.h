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
#include <algorithm>
#include <cctype>
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
    optimizer_type = NormalizeOptimizerType_(this->ParseStrOr("OptimizerType", "StochasticReconfiguration"));
    max_iterations = static_cast<size_t>(this->ParseIntOr("MaxIterations", 10));
    learning_rate = this->ParseDoubleOr("LearningRate", 0.01);
    
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
    } else if (optimizer_type == "LBFGS") {
      lbfgs_history_size = ParseNonNegativeSizeTOr_("LBFGSHistorySize", 10);
      lbfgs_tolerance_grad = this->ParseDoubleOr("LBFGSToleranceGrad", 1e-5);
      lbfgs_tolerance_change = this->ParseDoubleOr("LBFGSToleranceChange", 1e-9);
      lbfgs_max_eval = ParseNonNegativeSizeTOr_("LBFGSMaxEval", 20);
      lbfgs_step_mode = ParseLBFGSStepMode_(this->ParseStrOr("LBFGSStepMode", "Fixed"));
      lbfgs_wolfe_c1 = this->ParseDoubleOr("LBFGSWolfeC1", 1e-4);
      lbfgs_wolfe_c2 = this->ParseDoubleOr("LBFGSWolfeC2", 0.9);
      lbfgs_min_step = this->ParseDoubleOr("LBFGSMinStep", 1e-8);
      lbfgs_max_step = this->ParseDoubleOr("LBFGSMaxStep", 1.0);
      lbfgs_min_curvature = this->ParseDoubleOr("LBFGSMinCurvature", 1e-12);
      lbfgs_use_damping = this->ParseBoolOr("LBFGSUseDamping", true);
      lbfgs_max_direction_norm = this->ParseDoubleOr("LBFGSMaxDirectionNorm", 1e3);
      lbfgs_allow_fallback_to_fixed_step = this->ParseBoolOr("LBFGSAllowFallbackToFixedStep", false);
      lbfgs_fallback_fixed_step_scale = this->ParseDoubleOr("LBFGSFallbackFixedStepScale", 0.2);
    }
    
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
    
    // Gradient clipping (optional)
    double tmp_val = 0.0;
    if (this->TryParseDouble("ClipNorm", tmp_val)) clip_norm = tmp_val; else clip_norm.reset();
    if (this->TryParseDouble("ClipValue", tmp_val)) clip_value = tmp_val; else clip_value.reset();

    // Step selector parameters (optional; default: disabled)
    initial_step_selector_enabled = this->ParseBoolOr("InitialStepSelectorEnabled", false);
    initial_step_selector_max_line_search_steps =
        ParseNonNegativeSizeTOr_("InitialStepSelectorMaxLineSearchSteps", 3);
    initial_step_selector_enable_in_deterministic =
        this->ParseBoolOr("InitialStepSelectorEnableInDeterministic", false);
    auto_step_selector_enabled = this->ParseBoolOr("AutoStepSelectorEnabled", false);
    auto_step_selector_every_n_steps =
        ParseNonNegativeSizeTOr_("AutoStepSelectorEveryNSteps", 10);
    auto_step_selector_phase_switch_ratio = this->ParseDoubleOr("AutoStepSelectorPhaseSwitchRatio", 0.3);
    auto_step_selector_enable_in_deterministic =
        this->ParseBoolOr("AutoStepSelectorEnableInDeterministic", false);

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

    // Checkpointing (optional; default: disabled)
    checkpoint_every_n_steps = static_cast<size_t>(this->ParseIntOr("CheckpointEveryNSteps", 0));
    checkpoint_base_path = this->ParseStrOr("CheckpointBasePath", "");

    ValidateOptimizerConfig_();

    // Parse IO configuration
    io_params.Parse(*this);
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
  heisenberg_params::IOParams io_params; 
  
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

  // L-BFGS parameters
  size_t lbfgs_history_size = 10;
  double lbfgs_tolerance_grad = 1e-5;
  double lbfgs_tolerance_change = 1e-9;
  size_t lbfgs_max_eval = 20;
  qlpeps::LBFGSStepMode lbfgs_step_mode = qlpeps::LBFGSStepMode::kFixed;
  double lbfgs_wolfe_c1 = 1e-4;
  double lbfgs_wolfe_c2 = 0.9;
  double lbfgs_min_step = 1e-8;
  double lbfgs_max_step = 1.0;
  double lbfgs_min_curvature = 1e-12;
  bool lbfgs_use_damping = true;
  double lbfgs_max_direction_norm = 1e3;
  bool lbfgs_allow_fallback_to_fixed_step = false;
  double lbfgs_fallback_fixed_step_scale = 0.2;
  
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

  // Step selectors
  bool initial_step_selector_enabled = false;
  size_t initial_step_selector_max_line_search_steps = 3;
  bool initial_step_selector_enable_in_deterministic = false;
  bool auto_step_selector_enabled = false;
  size_t auto_step_selector_every_n_steps = 10;
  double auto_step_selector_phase_switch_ratio = 0.3;
  bool auto_step_selector_enable_in_deterministic = false;

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

  // Checkpointing
  size_t checkpoint_every_n_steps = 0;
  std::string checkpoint_base_path;

  /**
   * @brief Create qlpeps::VMCPEPSOptimizerParams with modern optimizer support
   */
  qlpeps::VMCPEPSOptimizerParams CreateVMCOptimizerParams(int rank = 0) {
    // Create Monte Carlo parameters
    size_t N = physical_params.Lx * physical_params.Ly;
    size_t spin_up_sites = N / 2;
    qlpeps::OccupancyNum occupancy = {spin_up_sites, N - spin_up_sites};
    qlpeps::Configuration config(physical_params.Ly, physical_params.Lx, occupancy);
    bool warmed_up = false;
    if (!io_params.configuration_load_dir.empty()) {
      warmed_up = config.Load(io_params.configuration_load_dir, static_cast<size_t>(rank));
    }

    qlpeps::MonteCarloParams mc_params_obj(
      mc_params.MC_total_samples,
      mc_params.WarmUp,
      mc_params.MCLocalUpdateSweepsBetweenSample,
      config,
      warmed_up,
      io_params.configuration_dump_dir
    );

    const qlpeps::PEPSParams peps_params_obj = heisenberg_params::CreatePEPSParams(
        physical_params.BoundaryCondition, bmps_params, *this);

    qlpeps::OptimizerParams opt_params = CreateOptimizerParams();

    return qlpeps::VMCPEPSOptimizerParams(opt_params, mc_params_obj, peps_params_obj);
  }

private:
  size_t ParseNonNegativeSizeTOr_(const std::string &key, int default_value) {
    const int value = this->ParseIntOr(key, default_value);
    if (value < 0) {
      throw std::invalid_argument(key + " must be >= 0");
    }
    return static_cast<size_t>(value);
  }

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
    base_params.initial_step_selector = qlpeps::InitialStepSelectorParams{
      initial_step_selector_enabled,
      initial_step_selector_max_line_search_steps,
      initial_step_selector_enable_in_deterministic
    };
    base_params.auto_step_selector = qlpeps::AutoStepSelectorParams{
      auto_step_selector_enabled,
      auto_step_selector_every_n_steps,
      auto_step_selector_phase_switch_ratio,
      auto_step_selector_enable_in_deterministic
    };
    
    // Create algorithm-specific parameters
    qlpeps::CheckpointParams ckpt_params{
      checkpoint_every_n_steps,
      checkpoint_base_path
    };
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
    } else if (optimizer_type == "LBFGS") {
      qlpeps::LBFGSParams lbfgs_params(
        lbfgs_history_size,
        lbfgs_tolerance_grad,
        lbfgs_tolerance_change,
        lbfgs_max_eval,
        lbfgs_step_mode,
        lbfgs_wolfe_c1,
        lbfgs_wolfe_c2,
        lbfgs_min_step,
        lbfgs_max_step,
        lbfgs_min_curvature,
        lbfgs_use_damping,
        lbfgs_max_direction_norm,
        lbfgs_allow_fallback_to_fixed_step,
        lbfgs_fallback_fixed_step_scale
      );
      return qlpeps::OptimizerParams(base_params, lbfgs_params, ckpt_params, spike_params);
    } else {
      throw std::invalid_argument("Unknown optimizer type: " + optimizer_type);
    }
  }

  static std::string ToLowerAscii_(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
  }

  static std::string NormalizeOptimizerType_(const std::string &value) {
    const std::string key = ToLowerAscii_(heisenberg_params::TrimAsciiWhitespace(value));
    if (key == "sr" || key == "stochasticreconfiguration") return "StochasticReconfiguration";
    if (key == "sgd") return "SGD";
    if (key == "adam") return "Adam";
    if (key == "adagrad") return "AdaGrad";
    if (key == "lbfgs" || key == "l-bfgs") return "LBFGS";
    return heisenberg_params::TrimAsciiWhitespace(value);
  }

  static qlpeps::LBFGSStepMode ParseLBFGSStepMode_(const std::string &value) {
    const std::string key = ToLowerAscii_(heisenberg_params::TrimAsciiWhitespace(value));
    if (key == "fixed" || key == "kfixed") return qlpeps::LBFGSStepMode::kFixed;
    if (key == "strongwolfe" || key == "kstrongwolfe") return qlpeps::LBFGSStepMode::kStrongWolfe;
    throw std::invalid_argument("LBFGSStepMode must be one of: Fixed, StrongWolfe, kFixed, kStrongWolfe");
  }

  void ValidateOptimizerConfig_() const {
    if (optimizer_type != "SGD" && optimizer_type != "Adam" &&
        optimizer_type != "AdaGrad" && optimizer_type != "StochasticReconfiguration" &&
        optimizer_type != "LBFGS") {
      throw std::invalid_argument("Unknown optimizer type: " + optimizer_type);
    }

    const bool any_selector_enabled = initial_step_selector_enabled || auto_step_selector_enabled;
    if (any_selector_enabled &&
        !(optimizer_type == "SGD" || optimizer_type == "StochasticReconfiguration")) {
      throw std::invalid_argument(
          "Step selectors only support SGD and StochasticReconfiguration");
    }
    if (any_selector_enabled && !lr_scheduler_type.empty()) {
      throw std::invalid_argument(
          "Step selectors cannot be used together with LRScheduler");
    }
    if (any_selector_enabled && learning_rate <= 0.0) {
      throw std::invalid_argument("Step selectors require a positive LearningRate");
    }
    if (initial_step_selector_enabled && initial_step_selector_max_line_search_steps == 0) {
      throw std::invalid_argument(
          "InitialStepSelectorMaxLineSearchSteps must be > 0 when InitialStepSelectorEnabled=true");
    }
    if (auto_step_selector_enabled) {
      if (auto_step_selector_every_n_steps == 0) {
        throw std::invalid_argument(
            "AutoStepSelectorEveryNSteps must be > 0 when AutoStepSelectorEnabled=true");
      }
      if (auto_step_selector_phase_switch_ratio < 0.0 ||
          auto_step_selector_phase_switch_ratio > 1.0) {
        throw std::invalid_argument(
            "AutoStepSelectorPhaseSwitchRatio must be within [0, 1]");
      }
    }

    if (optimizer_type != "LBFGS") {
      return;
    }
    if (lbfgs_history_size == 0) {
      throw std::invalid_argument("LBFGSHistorySize must be > 0");
    }
    if (lbfgs_step_mode == qlpeps::LBFGSStepMode::kStrongWolfe) {
      if (!(lbfgs_wolfe_c1 > 0.0 && lbfgs_wolfe_c1 < lbfgs_wolfe_c2 &&
            lbfgs_wolfe_c2 < 1.0)) {
        throw std::invalid_argument(
            "Strong-Wolfe constants must satisfy 0 < LBFGSWolfeC1 < LBFGSWolfeC2 < 1");
      }
      if (lbfgs_max_eval == 0) {
        throw std::invalid_argument("LBFGSMaxEval must be > 0 in StrongWolfe mode");
      }
      if (lbfgs_tolerance_grad < 0.0) {
        throw std::invalid_argument("LBFGSToleranceGrad must be >= 0 in StrongWolfe mode");
      }
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
