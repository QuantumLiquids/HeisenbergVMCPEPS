#include "qlten/qlten.h"
#include <ctime>
#include "qldouble.h"
#include "params_parser.h"
#include "myutil.h"
#include "dmrg_my_measure.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

using Link = std::pair<size_t, size_t>;

std::vector<Link> GenerateOBCTriangularNNLink(const size_t Lx, const size_t Ly) {
  using std::make_pair;
  size_t N = Lx * Ly;
  std::vector<Link> res;
  res.reserve(3 * N);

  for (size_t x = 0; x < Lx; x++) {
    for (size_t y = 0; y < Ly; y++) {
      size_t site0 = Ly * x + y;
      if (y < Ly - 1) {
        res.push_back(make_pair(site0, site0 + 1)); //vertical bond
      } else {
//        res.push_back(make_pair(site0 - (Ly - 1), site0)); //vertical winding bond, for cylinder
      }

      if (x < Lx - 1) {
        res.push_back(make_pair(site0, site0 + Ly)); //horizontal bond
        if (y < Ly - 1) {
          res.push_back(make_pair(site0 + 1, site0 + Ly)); //diagonal bond
        } else {
//          res.push_back(make_pair(site0 - (Ly - 1), site0 + Ly)); //diagonal bond, for cylinder
        }
      }
    }
  }
  return res;
}

int main(int argc, char *argv[]) {
//  namespace mpi = boost::mpi;
//  mpi::environment env;
//  mpi::communicator world;

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

  Timer one_site_timer("measure  one site operators");
  MeasureOneSiteOp(mps, kMpsPath, sz, "sz");
  cout << "measured one point function.<====" << endl;
  one_site_timer.PrintElapsed();

  std::vector<Link> links;
  links = GenerateOBCTriangularNNLink(Lx, Ly);

  std::vector<std::vector<size_t>> sites_set(links.size());
  for (size_t i = 0; i < links.size(); i++) {
    sites_set[i] = {links[i].first, links[i].second};
  }

  mps.Load();
  MeasureTwoSiteOp(mps, {sz, sz}, id, sites_set, "bondzz");
  MeasureTwoSiteOp(mps, {sp, sm}, id, sites_set, "bondpm");
  MeasureTwoSiteOp(mps, {sm, sp}, id, sites_set, "bondmp");

  std::vector<std::vector<size_t>> corr_sites_set(Lx / 2);
  for (size_t i = 0; i < corr_sites_set.size(); i++) {
    size_t site1 = Lx / 4 * Ly + Ly / 2;
    size_t site2 = (Lx / 4 + i + 1) * Ly + Ly / 2;
    corr_sites_set[i] = {site1, site2};
  }
  MeasureTwoSiteOp(mps, {sz, sz}, id, sites_set, "corrzz");
  MeasureTwoSiteOp(mps, {sp, sm}, id, sites_set, "corrpm");
  MeasureTwoSiteOp(mps, {sm, sp}, id, sites_set, "corrmp");

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
