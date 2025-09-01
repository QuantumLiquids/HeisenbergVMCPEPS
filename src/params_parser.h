//
// Created by haoxinwang on 28/09/2023.
//

#ifndef HEISENBERGVMCPEPS_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer_params.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_peps_params.h"
#include "common_params.h"

struct SimpleUpdateParams : public qlmps::CaseParamsParserBasic {
  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::NumericalParams numerical_params;
  double Tau;    ///< Time step for imaginary time evolution
  size_t Step;   ///< Number of simple update steps

  SimpleUpdateParams(const char *physics_file, const char *algorithm_file) : 
      qlmps::CaseParamsParserBasic(algorithm_file),
      physical_params(physics_file),
      numerical_params(algorithm_file) {
    Tau = ParseDouble("Tau");
    Step = ParseInt("Step");
  }
};

struct VMCUpdateParams : public qlmps::CaseParamsParserBasic {
  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  
  // VMC-specific parameters
  size_t CGMaxIter;                                    ///< Conjugate gradient max iterations
  double CGTol;                                        ///< Conjugate gradient tolerance
  int CGResidueRestart;                               ///< CG residue restart parameter
  double CGDiagShift;                                 ///< CG diagonal shift
  bool ReplicaTest;                                   ///< Whether to perform replica test
  int update_scheme;   ///< Wave function update scheme (legacy, now handled by OptimizerParams factory)
  std::vector<double> step_len;                       ///< Step lengths for optimization

  VMCUpdateParams(const char *physics_file, const char *algorithm_file) : 
      qlmps::CaseParamsParserBasic(algorithm_file),
      physical_params(physics_file),
      mc_params(algorithm_file),
      bmps_params(algorithm_file) {
    CGMaxIter = ParseInt("CGMaxIter");
    CGTol = ParseDouble("CGTol");
    CGResidueRestart = ParseInt("CGResidueRestart");
    CGDiagShift = ParseDouble("CGDiagShift");
    ReplicaTest = ParseBool("ReplicaTest");
    
    size_t update_times = ParseInt("UpdateNum");
    step_len = std::vector<double>(update_times);
    if (update_times > 0) {
      step_len[0] = ParseDouble("StepLengthFirst");
      double step_len_change = ParseDouble("StepLengthDecrease");
      for (size_t i = 1; i < update_times; i++) {
        step_len[i] = step_len[0] - i * step_len_change;
      }
    }
    update_scheme = ParseInt("UpdateScheme");
  }
  
  /**
   * @brief Create qlpeps::VMCPEPSOptimizerParams from these parameters
   */
  qlpeps::VMCPEPSOptimizerParams CreateVMCOptimizerParams() const {
    size_t N = physical_params.Lx * physical_params.Ly;
    
    // Create initial configuration with half spin up, half spin down
    qlpeps::Configuration config(physical_params.Ly, physical_params.Lx);
    size_t spin_up_sites = N / 2;
    for (size_t row = 0; row < physical_params.Ly; row++) {
      for (size_t col = 0; col < physical_params.Lx; col++) {
        size_t site_idx = row * physical_params.Lx + col;
        config({row, col}) = (site_idx < spin_up_sites) ? 1 : 0;
      }
    }
    
    qlpeps::MonteCarloParams mc_params_obj(
      mc_params.MC_samples,
      mc_params.WarmUp,
      mc_params.MCLocalUpdateSweepsBetweenSample,
      config,
      false  // not warmed up
    );
    
    qlpeps::PEPSParams peps_params_obj(
      bmps_params.CreateTruncatePara()
    );
    
    // Create optimizer parameters using factory method for Stochastic Reconfiguration
    qlpeps::OptimizerParams opt_params = qlpeps::OptimizerFactory::CreateStochasticReconfiguration(
      step_len.size(),  // max_iterations
      qlpeps::ConjugateGradientParams(CGMaxIter, CGTol, CGResidueRestart, CGDiagShift),
      !step_len.empty() ? step_len[0] : 0.01  // learning_rate
    );
    
    return qlpeps::VMCPEPSOptimizerParams(opt_params, mc_params_obj, peps_params_obj);
  }
};

/**
 * @brief Parameters for Monte Carlo measurement
 */
struct MCMeasurementParams : public qlmps::CaseParamsParserBasic {
  heisenberg_params::PhysicalParams physical_params;
  heisenberg_params::MonteCarloNumericalParams mc_params;
  heisenberg_params::BMPSParams bmps_params;
  
  MCMeasurementParams(const char *physics_file, const char *algorithm_file) : 
      qlmps::CaseParamsParserBasic(algorithm_file),
      physical_params(physics_file),
      mc_params(algorithm_file),
      bmps_params(algorithm_file) {
    // Measurement-specific parameters can be added here
    // For now, we use the common parameters
  }
  
  /**
   * @brief Create qlpeps::MCMeasurementParams from these parameters
   */
  qlpeps::MCMeasurementParams CreateMCMeasurementParams() const {
    size_t N = physical_params.Lx * physical_params.Ly;
    
    // Create initial configuration with half spin up, half spin down
    qlpeps::Configuration config(physical_params.Ly, physical_params.Lx);
    size_t spin_up_sites = N / 2;
    for (size_t row = 0; row < physical_params.Ly; row++) {
      for (size_t col = 0; col < physical_params.Lx; col++) {
        size_t site_idx = row * physical_params.Lx + col;
        config({row, col}) = (site_idx < spin_up_sites) ? 1 : 0;
      }
    }
    
    qlpeps::MonteCarloParams mc_params_obj(
      mc_params.MC_samples,
      mc_params.WarmUp,
      mc_params.MCLocalUpdateSweepsBetweenSample,
      config,
      false  // not warmed up
    );
    
    qlpeps::PEPSParams peps_params_obj(
      bmps_params.CreateTruncatePara()
    );
    
    return qlpeps::MCMeasurementParams(mc_params_obj, peps_params_obj);
  }
};

struct DMRGCaseParams : public qlmps::CaseParamsParserBasic {
  DMRGCaseParams(const char *pf) : CaseParamsParserBasic(pf) {
//    Geometry = ParseStr("Geometry");
    Ly = ParseInt("Ly");
    Lx = ParseInt("Lx");
    RemoveCorner = ParseBool("RemoveCorner");
//    J1 = ParseDouble("J1");
    J2 = ParseDouble("J2");
//    J3 = ParseDouble("J3");
//    Dzz = ParseDouble("Dzz");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
//    tau = ParseDouble("tau");
//    M = ParseInt("M");
//    ConvergeTolerance = ParseDouble("ConvergeTolerance");
//    TaylorOrder = ParseInt("TaylorOrder");
//    TaylorErr = ParseDouble("TaylorErr");
    Threads = ParseInt("Threads");
//    Perturbation=ParseDouble("Perturbation");
//    wavelength = ParseInt("wavelength");
    noise = ParseDoubleVec("noise");
//    SymmetryMode = ParseInt("SymmetryMode");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  bool RemoveCorner;
  double J1;
  double J2;
  double J3;
  double Dzz;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  double tau;
  size_t M;
  double ConvergeTolerance;
  size_t TaylorOrder;
  double TaylorErr;
  size_t Threads;
  double Perturbation;
  size_t wavelength;
  std::vector<double> noise;
  size_t SymmetryMode;//useless upto now
};

#endif //HEISENBERGVMCPEPS_PARAMS_PARSER_H
