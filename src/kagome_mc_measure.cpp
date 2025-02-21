//
// Created by haoxinwang on 04/11/2023.
//

#include "qlpeps/two_dim_tn/peps/square_lattice_peps.h"
#include "qlpeps/algorithm/vmc_update/monte_carlo_measurement.h"
#include "spin_onehalf_heisenberg_kagome_model_sqrpeps_solver.h"
#include "./qldouble.h"
#include "./params_parser.h"
#include "kagome_hei_model_combined_tps_sample.h"

using namespace qlpeps;

using MCUpdater = KagomeCombinedTPSSampleLoaclFlip<TenElemT, U1QN>;

Configuration GenerateInitialConfigurationInSmoothBoundary(size_t ly, size_t lx) {
  size_t sys_ly = 2 * ly, sys_lx = 2 * lx;
  std::vector<std::vector<size_t>> activates(sys_ly, std::vector<size_t>(sys_lx, 0));
  size_t sz_int = 0;
  for (size_t y = 0; y < sys_ly - 1; y++) {
    for (size_t x = 0; x < sys_lx - 1; x++) {
      if ((x & 1) && (y & 1)) {
        activates[y][x] = 0;
      } else {
        activates[y][x] = sz_int % 2;
        sz_int += 1;
      }
    }
  }
  Configuration config(ly, lx);
  for (size_t row = 0; row < ly; row++) {
    for (size_t col = 0; col < lx; col++) {
      config({row, col}) = activates[2 * row][2 * col]
          + 2 * activates[2 * row + 1][2 * col]
          + 4 * activates[2 * row][2 * col + 1];
    }
  }
  return config;
}

std::vector<std::array<size_t, 4>> GenerateEnergyBondInfo(size_t ly, size_t lx) {
  std::vector<std::array<size_t, 4>> bonds;
  bonds.reserve(6 * lx * ly);
  std::array<size_t, 4> bond;
  for (size_t row = 0; row < ly - 1; row++) {
    for (size_t col = 0; col < lx - 1; col++) {
      //[y1,x1, y2, x2]
      //shear to rectangular
      bond = {2 * row, 2 * col, 2 * row + 1, 2 * col};
      bonds.push_back(bond);
      bond = {2 * row + 1, 2 * col, 2 * row, 2 * col + 1};
      bonds.push_back(bond);
      bond = {2 * row, 2 * col, 2 * row, 2 * col + 1};
      bonds.push_back(bond);
      bond = {2 * row, 2 * col + 1, 2 * row, 2 * col + 2};
      bonds.push_back(bond);
    }
    bond = {2 * row, 2 * (ly - 1), 2 * row + 1, 2 * (ly - 1)};
    bonds.push_back(bond);
    for (size_t col = 0; col < lx - 1; col++) {
      bond = {2 * row + 2, 2 * col + 1, 2 * row + 1, 2 * col + 2};
      bonds.push_back(bond);
    }
  }
  for (size_t col = 0; col < lx - 1; col++) {
    bond = {2 * (ly - 1), 2 * col, 2 * (ly - 1), 2 * col + 1};
    bonds.push_back(bond);
    bond = {2 * (ly - 1), 2 * col + 1, 2 * (ly - 1), 2 * col + 2};
    bonds.push_back(bond);
  }

  for (size_t col = 0; col < lx; col++) {
    for (size_t row = 0; row < ly - 1; row++) {
      bond = {2 * row + 1, 2 * col, 2 * row + 2, 2 * col};
      bonds.push_back(bond);
    }
  }
  return bonds;
}

void DumpBondInfo(size_t ly, size_t lx, std::string &basename) {
  std::vector<std::array<size_t, 4>> bonds = GenerateEnergyBondInfo(ly, lx);
  auto file = basename + ".json";
  std::ofstream ofs(file);
  ofs << "[\n";
  for (size_t i = 0; i < bonds.size(); i++) {
    std::array<size_t, 4> &bond = bonds[i];
    ofs << "[" << bond[0] << ", " << bond[1] << ", " << bond[2] << ", " << bond[3] << "]";
    if (i < bonds.size() - 1)
      ofs << ",\n";
    else
      ofs << "\n";
  }
  ofs << "]";
  ofs.close();
}

int main(int argc, char **argv) {
 MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  VMCUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  size_t N = params.Lx * params.Ly;
  std::vector<size_t> occupation_num(8, 0);

  size_t occupation = N / 6;  //no 0 and 7;
  size_t remain_occupation = N - 6 * occupation;
  for (size_t i = 1; i < 7; i++) {
    occupation_num[i] = occupation;
  }
  if (remain_occupation > 0) {
    occupation_num[1] += remain_occupation / 2;
    occupation_num[3] += remain_occupation / 2;
  }
  if (remain_occupation % 2 == 1) {
    std::cout << "lattice site number = odd, not support." << std::endl;
    exit(-1);
  }

  qlpeps::MCMeasurementPara measurement_para;
  if (params.RemoveCorner) {
    Configuration init_config = GenerateInitialConfigurationInSmoothBoundary(params.Ly, params.Lx);
    measurement_para = qlpeps::MCMeasurementPara(
        BMPSTruncatePara(params.Db_min, params.Db_max,
                         params.TruncErr,
                         params.MPSCompressScheme,
                         std::make_optional<double>(params.TruncErr),
                         std::make_optional<size_t>(10)),
        params.MC_samples, params.WarmUp,
        params.MCLocalUpdateSweepsBetweenSample,
        init_config);
  } else {
    measurement_para = qlpeps::MCMeasurementPara(
        BMPSTruncatePara(params.Db_min, params.Db_max,
                         params.TruncErr,
                         params.MPSCompressScheme,
                         std::make_optional<double>(params.TruncErr),
                         std::make_optional<size_t>(10)),
        params.MC_samples, params.WarmUp,
        params.MCLocalUpdateSweepsBetweenSample,
        occupation_num,
        params.Ly, params.Lx);
  }

  using Model = KagomeSpinOneHalfHeisenbergOnSquarePEPSSolver<TenElemT, U1QN>;
  MonteCarloMeasurementExecutor<TenElemT, U1QN, MCUpdater, Model> *executor(nullptr);
  if (qlmps::IsPathExist(measurement_para.wavefunction_path)) {
    executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, MCUpdater, Model>(measurement_para,
                                                                                    params.Ly, params.Lx,
                                                                                    world);
  } else {  //after simple update, load from PEPS tensors
    SquareLatticePEPS<QLTEN_Double, U1QN> peps(pb_out, 2 * params.Ly, 2 * params.Lx);
    if (!peps.Load(peps_path)) {
      std::cout << "Loading simple updated PEPS files is broken." << std::endl;
      exit(-2);
    };
    SplitIndexTPS<QLTEN_Double, U1QN> split_idx_tps = KagomeSquarePEPSToSplitIndexTPS(peps);
    if (!split_idx_tps.IsBondDimensionEven()) {
      std::cout << "Warning: Split Index TPS bond dimension  is not even!" << std::endl;
    }
    executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, MCUpdater, Model>(measurement_para,
                                                                                    split_idx_tps,
                                                                                    world);
  }

  if (params.ReplicaTest) {
//    executor->ReplicaTest();
  } else {
    executor->Execute();
  }

  delete executor;
  std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);
  DumpBondInfo(params.Ly, params.Lx, bondinfo_filename);
  return 0;
}
