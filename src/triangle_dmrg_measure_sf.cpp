// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-02-10
*
* Description: Measurement of spin structure factor on the triangle heisenberg models.
*/

#include "qlten/qlten.h"
#include <ctime>
#include "qldouble.h"
#include "params_parser.h"
#include "myutil.h"
#include "dmrg_my_measure.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  DMRGCaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = Lx * Ly;
  clock_t startTime, endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  Tensor sz = Tensor({pb_in, pb_out});
  Tensor sp = Tensor({pb_in, pb_out});
  Tensor sm = Tensor({pb_in, pb_out});
  Tensor id = Tensor({pb_in, pb_out});
  sz({0, 0}) = 0.5;
  sz({1, 1}) = -0.5;
  sp({0, 1}) = 1.0;
  sm({1, 0}) = 1.0;
  id({0, 0}) = 1.0;
  id({1, 1}) = 1.0;
  const SiteVec<TenElemT, U1QN> sites = SiteVec<TenElemT, U1QN>(N, pb_out);

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  /******** define the measure_tasks ********/
  std::vector<MeasureGroupTask> measure_tasks;
  measure_tasks.reserve(N);
  size_t begin_x = 0;
  size_t end_x = Lx;
  for (size_t i = begin_x * Ly; i < end_x * Ly; i++) {
    const size_t site1 = i;
    std::vector<size_t> site2;
    site2.reserve(N - i);
    for (size_t j = i + 1; j < end_x * Ly; j++) {
      site2.push_back(j);
    }
    measure_tasks.push_back(MeasureGroupTask(site1, site2));
  }
  mps.Load();
  Timer measure_struct_factor_timer("measure structure factor");
  MeasureTwoSiteOp(mps, sz, sz, measure_tasks, "zzsf", world);
  MeasureTwoSiteOp(mps, sp, sm, measure_tasks, "pmsf", world);
  MeasureTwoSiteOp(mps, sm, sp, measure_tasks, "mpsf", world);
  measure_struct_factor_timer.PrintElapsed();
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
