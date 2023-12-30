//
// Created by haoxinwang on 28/09/2023.
//

#ifndef HEISENBERGVMCPEPS_PARAMS_PARSER_H
#define HEISENBERGVMCPEPS_PARAMS_PARSER_H

#include "gqmps2/case_params_parser.h"
#include "gqpeps/algorithm/vmc_update/vmc_peps.h"

struct SimpleUpdateParams : public gqmps2::CaseParamsParserBasic {
  SimpleUpdateParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    RemoveCorner = ParseBool("RemoveCorner");
    J2 = ParseDouble("J2");
    TruncErr = ParseDouble("TruncErr");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    Tau = ParseDouble("Tau");
    Step = ParseInt("Step");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Ly;
  size_t Lx;
  bool RemoveCorner;
  double J2;
  double TruncErr;
  size_t Dmin;
  size_t Dmax;
  double Tau;
  size_t Step;
  size_t ThreadNum;
};

struct VMCUpdateParams : public gqmps2::CaseParamsParserBasic {
  VMCUpdateParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    J2 = ParseDouble("J2");
    Db_min = ParseInt("Dbmps_min");
    Db_max = ParseInt("Dbmps_max");
    TruncErr = ParseDouble("TruncErr");
    MC_samples = ParseInt("MC_samples");
    WarmUp = ParseInt("WarmUp");
    MCLocalUpdateSweepsBetweenSample = ParseInt("MCLocalUpdateSweepsBetweenSample");
    CGMaxIter = ParseInt("CGMaxIter");
    CGTol = ParseDouble("CGTol");
    CGResidueRestart = ParseInt("CGResidueRestart");
    CGDiagShift = ParseDouble("CGDiagShift");
    ReplicaTest = ParseBool("ReplicaTest");
    MPSCompressScheme = static_cast<gqpeps::CompressMPSScheme>(ParseInt("MPSCompressScheme"));
    RemoveCorner = ParseBool("RemoveCorner");
    size_t update_times = ParseInt("UpdateNum");
    step_len = std::vector<double>(update_times);
    if (update_times > 0) {
      step_len[0] = ParseDouble("StepLengthFirst");
      double step_len_change = ParseDouble("StepLengthDecrease");
      for (size_t i = 1; i < update_times; i++) {
        step_len[i] = step_len[0] - i * step_len_change;
      }
    }
    update_scheme = (gqpeps::WAVEFUNCTION_UPDATE_SCHEME) ParseInt("UpdateScheme");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Ly;
  size_t Lx;
  double J2;
  size_t Db_min;
  size_t Db_max;
  size_t TruncErr;
  size_t MC_samples;
  size_t WarmUp;
  size_t MCLocalUpdateSweepsBetweenSample;
  size_t CGMaxIter;
  double CGTol;
  int CGResidueRestart;
  double CGDiagShift;
  bool ReplicaTest;
  gqpeps::CompressMPSScheme MPSCompressScheme;
  bool RemoveCorner;
  gqpeps::WAVEFUNCTION_UPDATE_SCHEME update_scheme;
  std::vector<double> step_len;
  size_t ThreadNum;
};

struct DMRGCaseParams : public gqmps2::CaseParamsParserBasic {
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
