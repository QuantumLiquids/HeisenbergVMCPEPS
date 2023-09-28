/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2023-09-28
*
* Description: Simple Update for Heisenberg model.
*/


//#define PLAIN_TRANSPOSE 1

#include "gqpeps/algorithm/simple_update/triangle_nn_on_sqr_peps_simple_update.h"
#include "gqmps2/case_params_parser.h"
#include "./gqdouble.h"

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

std::string peps_path = "peps";


int main(int argc, char **argv) {
  SimpleUpdateParams params(argv[1]);

  Tensor did = Tensor({pb_in, pb_out});
  Tensor dsz = Tensor({pb_in, pb_out});
  Tensor dsp = Tensor({pb_in, pb_out});
  Tensor dsm = Tensor({pb_in, pb_out});
  Tensor ham_hei_nn = Tensor({pb_in, pb_out, pb_in, pb_out});

  did({0, 0}) = 1;
  did({1, 1}) = 1;
  dsz({0, 0}) = 0.5;
  dsz({1, 1}) = -0.5;
  dsp({0, 1}) = 1;
  dsm({1, 0}) = 1;

  ham_hei_nn({0, 0, 0, 0}) = 0.25;
  ham_hei_nn({1, 1, 1, 1}) = 0.25;
  ham_hei_nn({1, 1, 0, 0}) = -0.25;
  ham_hei_nn({0, 0, 1, 1}) = -0.25;
  ham_hei_nn({0, 1, 1, 0}) = 0.5;
  ham_hei_nn({1, 0, 0, 1}) = 0.5;

  Tensor ham_hei_tri_terms[3];
  for (size_t i = 0; i < 3; i++) {
    ham_hei_tri_terms[i] = Tensor({pb_in, pb_out, pb_in, pb_out, pb_in, pb_out});
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[0]({0, 0, 0, 0, i, i}) = 0.25;
    ham_hei_tri_terms[0]({1, 1, 1, 1, i, i}) = 0.25;
    ham_hei_tri_terms[0]({1, 1, 0, 0, i, i}) = -0.25;
    ham_hei_tri_terms[0]({0, 0, 1, 1, i, i}) = -0.25;
    ham_hei_tri_terms[0]({0, 1, 1, 0, i, i}) = 0.5;
    ham_hei_tri_terms[0]({1, 0, 0, 1, i, i}) = 0.5;
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[1]({0, 0, i, i, 0, 0}) = 0.25;
    ham_hei_tri_terms[1]({1, 1, i, i, 1, 1}) = 0.25;
    ham_hei_tri_terms[1]({1, 1, i, i, 0, 0}) = -0.25;
    ham_hei_tri_terms[1]({0, 0, i, i, 1, 1}) = -0.25;
    ham_hei_tri_terms[1]({0, 1, i, i, 1, 0}) = 0.5;
    ham_hei_tri_terms[1]({1, 0, i, i, 0, 1}) = 0.5;
  }

  for (size_t i = 0; i < 2; i++) {
    ham_hei_tri_terms[2]({i, i, 0, 0, 0, 0}) = 0.25;
    ham_hei_tri_terms[2]({i, i, 1, 1, 1, 1}) = 0.25;
    ham_hei_tri_terms[2]({i, i, 1, 1, 0, 0}) = -0.25;
    ham_hei_tri_terms[2]({i, i, 0, 0, 1, 1}) = -0.25;
    ham_hei_tri_terms[2]({i, i, 0, 1, 1, 0}) = 0.5;
    ham_hei_tri_terms[2]({i, i, 1, 0, 0, 1}) = 0.5;
  }
  Tensor ham_hei_tri = ham_hei_tri_terms[0] + ham_hei_tri_terms[1] + ham_hei_tri_terms[2];

  gqten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.ThreadNum);

  gqpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);

  gqpeps::SquareLatticePEPS<TenElemT, U1QN> peps0(pb_out, params.Ly, params.Lx);
  if (gqmps2::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    std::vector<std::vector<size_t>> activates(params.Ly, std::vector<size_t>(params.Lx));
    for (size_t y = 0; y < params.Ly; y++) {
      for (size_t x = 0; x < params.Lx; x++) {
        size_t sz_int = x + y;
        activates[y][x] = sz_int % 2;
      }
    }
    peps0.Initial(activates);
  }
  auto su_exe = new gqpeps::TriangleNNModelSquarePEPSSimpleUpdateExecutor(update_para, peps0,
                                                                          ham_hei_nn,
                                                                          ham_hei_tri);
  su_exe->Execute();
  auto tps = gqpeps::TPS<TenElemT, U1QN>(su_exe->GetPEPS());
  tps.Dump();
  su_exe->DumpResult(peps_path, true);
  return 0;
}