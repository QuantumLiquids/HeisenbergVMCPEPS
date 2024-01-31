/**
 * 2 processes parallel
 */
#include "qlmps/qlmps.h"
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

std::vector<Link> GenerateOBCKagomeNNLinkWithCorner(const size_t Lx, const size_t Ly) {
  using std::make_pair;
  size_t N = 3 * Lx * Ly;
  std::vector<Link> res;
  res.reserve(2 * N);
  for (size_t x = 0; x < Lx; x++) {
    for (size_t y = 0; y < Ly; y++) {
      size_t site0 = (3 * Ly) * x + 2 * y;
      size_t site1 = site0 + 1;
      size_t site2 = (3 * Ly) * x + 2 * Ly + y;
      res.push_back(make_pair(site0, site1));
      res.push_back(make_pair(site0, site2));
      res.push_back(make_pair(site1, site2));
      if (y < Ly - 1) {
        res.push_back(make_pair(site1, site1 + 1));
      }
      if (x < Lx - 1) {
        res.push_back(make_pair(site2, (3 * Ly) * (x + 1) + 2 * y));
      }
      if (y > 0 && x < Lx - 1) {
        res.push_back(make_pair(site2, (3 * Ly) * (x + 1) + 2 * y - 1));
      }
    }
  }
  return res;
}

std::vector<Link> GenerateOBCKagomeNNLinkSmoothBC(const size_t Lx, const size_t Ly) {
  using std::make_pair;
  size_t N = 3 * Lx * Ly - Lx - Ly;
  std::vector<Link> res;
  res.reserve(2 * N);
  for (size_t x = 0; x < Lx - 1; x++) {
    for (size_t y = 0; y < Ly - 1; y++) {
      size_t site0 = x * (3 * Ly - 1) + 2 * y;
      size_t site1 = site0 + 1;
      size_t site2 = x * (3 * Ly - 1) + (2 * Ly - 1) + y;
      res.push_back(make_pair(site0, site1));
      res.push_back(make_pair(site0, site2));
      res.push_back(make_pair(site1, site2));
      res.push_back(make_pair(site1, site1 + 1));
      res.push_back(make_pair(site2, (x + 1) * (3 * Ly - 1) + 2 * y));
      if (y > 0) {
        res.push_back(make_pair(site2, (x + 1) * (3 * Ly - 1) + 2 * y - 1));
      }
    }
  }

  //last row :
  for (size_t x = 0; x < Lx - 1; x++) {
    size_t y = Ly - 1;
    size_t site0 = x * (3 * Ly - 1) + 2 * y;
    size_t site2 = x * (3 * Ly - 1) + (2 * Ly - 1) + y;
    res.push_back(make_pair(site0, site2));
    res.push_back(make_pair(site2, (x + 1) * (3 * Ly - 1) + 2 * y));
    res.push_back(make_pair(site2, (x + 1) * (3 * Ly - 1) + 2 * y - 1));
  }
  //last column:
  for (size_t y = 0; y < Ly - 1; y++) {
    size_t x = Lx - 1;
    size_t site0 = x * (3 * Ly - 1) + 2 * y;
    size_t site1 = site0 + 1;
    res.push_back(make_pair(site0, site1));
    res.push_back(make_pair(site1, site1 + 1));
  }
  for (Link &link : res) {
    std::cout << "[ " << link.first << ", " << link.second << " ]" << std::endl;
  }
  return res;
}

int main(int argc, char *argv[]) {
//  namespace mpi = boost::mpi;
//  mpi::environment env;
//  mpi::communicator world;

  DMRGCaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N;
  if (params.RemoveCorner) {
    N = 3 * Lx * Ly - Lx - Ly;
  } else {
    N = 3 * Lx * Ly;
  }
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
  if (params.RemoveCorner) {
    links = GenerateOBCKagomeNNLinkSmoothBC(Lx, Ly);
  } else {
    links = GenerateOBCKagomeNNLinkWithCorner(Lx, Ly);
  }

  std::vector<std::vector<size_t>> sites_set(links.size());
  for (size_t i = 0; i < links.size(); i++) {
    sites_set[i] = {links[i].first, links[i].second};
  }

  mps.Load();
  MeasureTwoSiteOp(mps, {sz, sz}, id, sites_set, "bondzz");
  MeasureTwoSiteOp(mps, {sp, sm}, id, sites_set, "bondpm");
  MeasureTwoSiteOp(mps, {sm, sp}, id, sites_set, "bondmp");

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}