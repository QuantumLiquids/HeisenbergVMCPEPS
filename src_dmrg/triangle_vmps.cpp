//
// Created by haoxinwang on 28/12/2023.
//
#include <iostream>
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "../src/params_parser.h"
#include "../src/qldouble.h"
#include "../src/myutil.h"

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

std::vector<Link> GenerateOBCTriangularJ2Link(const size_t Lx, const size_t Ly) {
  using std::make_pair;
  size_t N = Lx * Ly;
  std::vector<Link> res;
  res.reserve(3 * N);

  for (size_t x = 0; x < Lx - 1; x++) {
    for (size_t y = 0; y < Ly - 1; y++) {
      size_t site0 = Ly * x + y;
      size_t site1 = Ly * (x + 1) + (y + 1);
      res.push_back(make_pair(site0, site1));
    }
  }
  for (size_t x = 0; x < Lx - 2; x++) {
    for (size_t y = 0; y < Ly - 1; y++) {
      size_t site0 = Ly * x + (y + 1);
      size_t site1 = Ly * (x + 2) + y;
      res.push_back(make_pair(site0, site1));
    }
  }
  for (size_t x = 0; x < Lx - 1; x++) {
    for (size_t y = 0; y < Ly - 2; y++) {
      size_t site0 = Ly * (x + 1) + y;
      size_t site1 = Ly * x + (y + 2);
      res.push_back(make_pair(site0, site1));
    }
  }
  return res;
}

int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);

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

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  double e0(0.0); //energy

  const SiteVec<TenElemT, QNT> sites(N, pb_out);
  qlmps::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);

  Tensor sz = Tensor({pb_in, pb_out});
  Tensor sp = Tensor({pb_in, pb_out});
  Tensor sm = Tensor({pb_in, pb_out});
  sz({0, 0}) = 0.5;
  sz({1, 1}) = -0.5;
  sp({0, 1}) = 1.0;
  sm({1, 0}) = 1.0;
  std::vector<Link> links = GenerateOBCTriangularNNLink(Lx, Ly);
  std::vector<Link> links_j2 = GenerateOBCTriangularJ2Link(Lx, Ly);

  for (auto &link : links) {
    mpo_gen.AddTerm(1.0, sz, link.first, sz, link.second);
    mpo_gen.AddTerm(0.5, sp, link.first, sm, link.second);
    mpo_gen.AddTerm(0.5, sm, link.first, sp, link.second);
  }
  for (auto &link : links_j2) {
    mpo_gen.AddTerm(1.0 * params.J2, sz, link.first, sz, link.second);
    mpo_gen.AddTerm(0.5 * params.J2, sp, link.first, sm, link.second);
    mpo_gen.AddTerm(0.5 * params.J2, sm, link.first, sp, link.second);
  }
  auto mpo = mpo_gen.Gen();

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N);

  for (size_t i = 0; i < N; ++i) {
    stat_labs[i] = i % 2;
  }

  if (rank == 0) {
    if (IsPathExist(kMpsPath)) {
      if (N == GetNumofMps()) {
        cout << "The number of mps files is consistent with mps size." << endl;
        cout << "Directly use mps from files." << endl;
      } else {
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(sweep_params.mps_path, true);
      }
    } else {
      qlmps::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(sweep_params.mps_path, true);
    }
  }

  if (!has_bond_dimension_parameter) {
    e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, comm);
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
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          D, D, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i]),
          params.noise
      );
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, comm);
    }
  }
  std::cout << "E0/site: " << e0 / N << std::endl;

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  MPI_Finalize();
  return 0;
}
