/**
 * 2 processes parallel
 */
#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <ctime>
#include "gqdouble.h"
#include "params_parser.h"
#include "myutil.h"
#include "dmrg_my_measure.h"


using namespace gqmps2;
using namespace gqten;
using namespace std;

int main(int argc, char *argv[]) {
//  namespace mpi = boost::mpi;
//  mpi::environment env;
//  mpi::communicator world;

  DMRGCaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = 3 * Lx * params.Ly;

  clock_t startTime, endTime;
  startTime = clock();

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);


  Tensor sz = Tensor({pb_in, pb_out});
  Tensor sp = Tensor({pb_in, pb_out});
  Tensor sm = Tensor({pb_in, pb_out});
  sz({0, 0}) = 0.5;
  sz({1, 1}) = -0.5;
  sp({0, 1}) = 1.0;
  sm({1, 0}) = 1.0;
  const SiteVec<TenElemT, U1QN> sites = SiteVec<TenElemT, U1QN>(N, pb_out);


  using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  Timer one_site_timer("measure  one site operators");
  MeasureOneSiteOp(mps, kMpsPath, sz, "sz");
  cout << "measured one point function.<====" << endl;
  one_site_timer.PrintElapsed();


  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}