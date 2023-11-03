//
// Created by haoxinwang on 02/11/2023.
//

#ifndef HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
#define HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H

#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "spin_onehalf_heisenberg_kagome_measurement_solver.h"

namespace gqpeps {
using namespace gqten;


template<typename TenElemT, typename QNT, typename MeasurementSolver>
class KagomeMeasurementExecutor : public Executor {
 public:
  using Tensor = GQTensor<TenElemT, QNT>;
  using TPST = TPS<TenElemT, QNT>;
  using SITPST = SplitIndexTPS<TenElemT, QNT>;
  using IndexT = Index<QNT>;

  //Load Data from path
  KagomeMeasurementExecutor(const VMCOptimizePara &optimize_para,
                            const size_t ly, const size_t lx,
                            const boost::mpi::communicator &world,
                            const MeasurementSolver &solver = MeasurementSolver());

  void Execute(void) override;

  void LoadTenData(void);

  void LoadTenData(const std::string &tps_path);

  void DumpData();

  void DumpData(const std::string &tps_path);


  VMCOptimizePara optimize_para;

 private:
  void ReserveSamplesDataSpace_();

  void PrintExecutorInfo_(void);

  void Measure_(void);

  size_t MCSweep_(void);

  void WarmUp_(void);

  void MeasureSample_(void);

  void GatherStatistic_(void);

  boost::mpi::communicator world_;

  size_t lx_; //cols
  size_t ly_; //rows

  SITPST split_index_tps_;

  TPSSample<TenElemT, QNT> tps_sample_;

  std::uniform_real_distribution<double> u_double_;

  bool warm_up_;

  MeasurementSolver measurement_solver_;

  // observable
  std::vector<TenElemT> energy_samples_;
  std::vector<size_t> center_configs_; //used to analyze the auto correlation.
  std::vector<std::vector<bool>> local_sz_samples_; // outside is the sample index, inner side is the lattice index.
  // the lattice site number = Lx * Ly * 3,  first the unit cell, then column idx, then row index.


};//KagomeMeasurementExecutor



template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::ReserveSamplesDataSpace_(void) {
  energy_samples_.reserve(optimize_para.mc_samples);
  center_configs_.reserve(optimize_para.mc_samples);
  local_sz_samples_.reserve(optimize_para.mc_samples);
  // add reserve if more measurements
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::MeasureSample_() {
  TenElemT energy;
  std::vector<bool> local_sz;
  energy = measurement_solver_(&split_index_tps_, &tps_sample_, local_sz);
  energy_samples_.push_back(energy);
  local_sz_samples_.push_back(local_sz);
  //add more measurement here and the definition of measurement solver
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::GatherStatistic_() {


}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void
KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::DumpData(const std::string &tps_path) {
  using gqmps2::IsPathExist;
  using gqmps2::CreatPath;
  tps_sample_.config.Dump(tps_path, world_.rank());
  std::string energy_raw_path = "energy_raw_data/";
  if (world_.rank() == kMasterProc)
    if (!IsPathExist(energy_raw_path)) {
      CreatPath(energy_raw_path);
    }
  world_.barrier();
  DumpVecData(energy_raw_path + "/energy" + std::to_string(world_.rank()), energy_samples_);

}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::KagomeMeasurementExecutor(
    const VMCOptimizePara &optimize_para,
    const size_t ly, const size_t lx,
    const boost::mpi::communicator &world,
    const MeasurementSolver &solver):
    world_(world), optimize_para(optimize_para), lx_(lx), ly_(ly),
    split_index_tps_(ly, lx), tps_sample_(ly, lx),
    u_double_(0, 1),
    energy_solver_(solver), warm_up_(false) {
  TPSSample<TenElemT, QNT>::trun_para = TruncatePara(optimize_para);
  random_engine.seed((size_t)
                         std::time(nullptr) + 10086 * world.rank());
  LoadTenData();
  ReserveSamplesDataSpace_();
  PrintExecutorInfo_();
  this->SetStatus(ExecutorStatus::INITED);
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::PrintExecutorInfo_(void) {
  if (world_.rank() == kMasterProc) {
    std::cout << std::left;  // Set left alignment for the output
    std::cout << "\n";
    std::cout << "=====> VARIATIONAL MONTE-CARLO PROGRAM FOR PEPS <=====" << "\n";
    std::cout << std::setw(30) << "System size (lx, ly):" << "(" << lx_ << ", " << ly_ << ")\n";
    std::cout << std::setw(30) << "PEPS bond dimension:" << split_index_tps_.GetMaxBondDimension() << "\n";
    std::cout << std::setw(30) << "BMPS bond dimension:" << optimize_para.bmps_trunc_para.D_min << "/"
              << optimize_para.bmps_trunc_para.D_max << "\n";
    std::cout << std::setw(30) << "Sampling numbers:" << optimize_para.mc_samples << "\n";

    std::cout << "=====> TECHNICAL PARAMETERS <=====" << "\n";
    std::cout << std::setw(40) << "The number of processors (including master):" << world_.size() << "\n";
    std::cout << std::setw(40) << "The number of threads per processor:" << hp_numeric::GetTensorManipulationThreads()
              << "\n";
  }
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::Execute(void) {
  SetStatus(ExecutorStatus::EXEING);
  Measure_();
  DumpData();
  SetStatus(ExecutorStatus::FINISH);
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::WarmUp_(void) {
  if (!warm_up_) {
    Timer warm_up_timer("warm_up");
    for (size_t sweep = 0; sweep < optimize_para.mc_warm_up_sweeps; sweep++) {
      MCSweep_();
    }
    double elasp_time = warm_up_timer.Elapsed();
    std::cout << "Proc " << std::setw(4) << world_.rank() << " warm-up completes T = " << elasp_time << "s."
              << std::endl;
    warm_up_ = true;
  }
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::LoadTenData(void) {
  LoadTenData(optimize_para.wavefunction_path);
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::LoadTenData(const std::string &tps_path) {
  if (!split_index_tps_.Load(tps_path)) {
    std::cout << "Loading TPS files fails." << std::endl;
    exit(-1);
  }
  Configuration config(ly_, lx_);
  bool load_config = config.Load(tps_path, world_.rank());
  if (load_config) {
    tps_sample_ = TPSSample<TenElemT, QNT>(split_index_tps_, config);
  } else {
    std::cout << "Loading configuration in rank " << world_.rank()
              << " fails. Random generate it and warm up."
              << std::endl;
    tps_sample_.RandomInit(split_index_tps_, optimize_para.occupancy_num, 10089 * world_.rank() + std::time(nullptr));
    WarmUp_();
  }
  warm_up_ = true;
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::DumpData(void) {
  DumpData(optimize_para.wavefunction_path);
}


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::Measure_(void) {
  size_t accept_num = 0;
  size_t bond_num = lx_ * (ly_ - 1) + ly_ * (lx_ - 1);
  for (size_t sweep = 0; sweep < optimize_para.mc_samples; sweep++) {
    accept_num += MCSweep_();
    MeasureSample_();
  }
  double accept_rate = double(accept_num) / double(bond_num * optimize_para.mc_samples);
  GatherStatistic_();
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
size_t KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::MCSweep_(void) {
  size_t flip_times;
  for (size_t i = 0; i < optimize_para.mc_sweeps_between_sample; i++) {
    flip_times = tps_sample_.MCCompressedKagomeLatticeLocalUpdateSweep(split_index_tps_, u_double_);
  }
  return flip_times;
}
}//gqpeps

#endif //HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
