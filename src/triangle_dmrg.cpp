//
// Created by haoxinwang on 28/12/2023.
//
#include <iostream>
#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "params_parser.h"
#include "gqdouble.h"
#include "myutil.h"

using namespace gqmps2;
using namespace gqten;
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
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  DMRGCaseParams params(argv[1]);
  const size_t Lx = params.Lx;
  const size_t Ly = params.Ly;
  size_t N = Lx * Ly;

  cout << "The total number of sites: " << N << endl;
  clock_t startTime, endTime;
  startTime = clock();

  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  gqmps2::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  double e0(0.0); //energy

  const SiteVec<TenElemT, U1QN> sites(N, pb_out);
  gqmps2::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);

  Tensor sz = Tensor({pb_in, pb_out});
  Tensor sp = Tensor({pb_in, pb_out});
  Tensor sm = Tensor({pb_in, pb_out});
  sz({0, 0}) = 0.5;
  sz({1, 1}) = -0.5;
  sp({0, 1}) = 1.0;
  sm({1, 0}) = 1.0;
  std::vector<Link> links = GenerateOBCTriangularNNLink(Lx, Ly);

  for (auto &link : links) {
    mpo_gen.AddTerm(1.0, sz, link.first, sz, link.second);
    mpo_gen.AddTerm(0.5, sp, link.first, sm, link.second);
    mpo_gen.AddTerm(0.5, sm, link.first, sp, link.second);
  }
  auto mpo = mpo_gen.Gen();

  using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N);

  for (size_t i = 0; i < N; ++i) {
    stat_labs[i] = i % 2;
  }

  if (world.rank() == 0) {
    if (IsPathExist(kMpsPath)) {
      if (N == GetNumofMps()) {
        cout << "The number of mps files is consistent with mps size." << endl;
        cout << "Directly use mps from files." << endl;
      } else {
        gqmps2::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(sweep_params.mps_path, true);
      }
    } else {
      gqmps2::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(sweep_params.mps_path, true);
    }
  }

  if (!has_bond_dimension_parameter) {
    e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
  } else {
    size_t DMRG_time = input_D_set.size();
    std::vector<size_t> MaxLanczIterSet(DMRG_time);
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        std::cout << MaxLanczIterSet[i] << ", ";
      }
      std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }

    for (size_t i = 0; i < DMRG_time; i++) {
      size_t D = input_D_set[i];
      std::cout << "D_max = " << D << std::endl;
      gqmps2::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          D, D, params.CutOff,
          gqmps2::LanczosParams(params.LanczErr, MaxLanczIterSet[i]),
          params.noise
      );
      e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
    }
  }
  std::cout << "E0/site: " << e0 / N << std::endl;

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
