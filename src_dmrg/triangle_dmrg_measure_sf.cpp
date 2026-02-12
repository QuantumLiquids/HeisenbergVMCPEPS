// SPDX-License-Identifier: LGPL-3.0-only

/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-02-10
*
* Description: Measurement of spin structure factor on the triangle heisenberg models.
*/

#include "qlten/qlten.h"
#include <ctime>
#include "../src/qldouble.h"
#include "dmrg_params.h"
#include "../src/myutil.h"
#include "dmrg_my_measure.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);

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
  const SiteVec<TenElemT, QNT> sites = SiteVec<TenElemT, QNT>(N, pb_out);

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
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
  MeasureTwoSiteOp(mps, sz, sz, measure_tasks, "zzsf", comm);
  MeasureTwoSiteOp(mps, sp, sm, measure_tasks, "pmsf", comm);
  MeasureTwoSiteOp(mps, sm, sp, measure_tasks, "mpsf", comm);
  measure_struct_factor_timer.PrintElapsed();
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  MPI_Finalize();
  return 0;
}
