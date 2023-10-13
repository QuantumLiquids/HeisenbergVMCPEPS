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
      Continue_from_VMC = ParseBool("Continue_from_VMC");
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
  bool Continue_from_VMC;
  gqpeps::WAVEFUNCTION_UPDATE_SCHEME update_scheme;
  std::vector<double> step_len;
  size_t ThreadNum;
};
#endif //HEISENBERGVMCPEPS_PARAMS_PARSER_H
