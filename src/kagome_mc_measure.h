//
// Created by haoxinwang on 02/11/2023.
//

#ifndef HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
#define HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H

#include "gqpeps/algorithm/vmc_update/vmc_peps.h"
#include "spin_onehalf_heisenberg_kagome_measurement_solver.h"

namespace gqpeps {
using namespace gqten;

std::vector<bool> KagomeConfig2Sz(
    const Configuration &config
) {
  std::vector<bool> local_sz;
  local_sz.reserve(config.size() * 3);
  for (size_t row = 0; row < config.rows(); row++) {
    for (size_t col = 0; col < config.cols(); col++) {
      size_t local_config = config({row, col});
      local_sz.push_back(local_config & 1); //left upper site
      local_sz.push_back(local_config >> 1 & 1);//lower site
      local_sz.push_back(local_config >> 2 & 1);//right site
    }
  }
  return local_sz;
}

///< sum (config1 * config2) in kagome lattice
size_t SpinConfigurationOverlap(
    const std::vector<bool> &sz1,
    const std::vector<bool> &sz2
) {
  size_t overlap(0);
  for (size_t i = 0; i < sz1.size(); i++) {
    overlap += sz1[i] && sz2[i];
  }
  return overlap;
}

///< 1/N * sum (sz1 * sz2)
double SpinConfigurationOverlap2(
    const std::vector<bool> &sz1,
    const std::vector<bool> &sz2
) {
  int overlap_sum(0);
  for (size_t i = 0; i < sz1.size(); i++) {
    overlap_sum += (2 * (int) sz1[i] - 1) * (2 * (int) sz2[i] - 1);
  }
  return double(overlap_sum) / sz1.size();
}

template<typename T>
std::vector<T> CalAutoCorrelation(
    const std::vector<T> &data,
    const T mean
) {
  const size_t res_len = 20; // I think enough long
  std::vector<T> res(res_len, T(0));
  for (size_t t = 0; t < res_len; t++) {
    T sum(0);
    for (size_t j = 0; j < data.size() - t; j++) {
      sum += data[j] * data[j + t];
    }
    res[t] = sum / (data.size() - t) - mean * mean;
  }
  return res;
}

template<typename T>
std::vector<T> AveListOfData(
    const std::vector<std::vector<T> > &data //outside idx: sample index; inside idx: something like site/bond
) {
  const size_t N = data[0].size();
  const size_t sample_size = data.size();
  std::vector<T> sum(N, T(0)), ave(N);
  for (size_t sample_idx = 0; sample_idx < sample_size; sample_idx++) {
    for (size_t i = 0; i < N; i++) {
      sum[i] += data[sample_idx][i];
    }
  }
  for (size_t i = 0; i < N; i++) {
    ave[i] = sum[i] / sample_size;
  }
  return ave;
}

std::vector<double> CalSpinAutoCorrelation(
    const std::vector<std::vector<bool>> &local_sz_samples
) {
  const size_t res_len = 20;
  const size_t N = local_sz_samples[0].size();// lattice size
  std::vector<double> res(res_len, 0.0);
  for (size_t t = 0; t < res_len; t++) {
    size_t overlap_sum(0);
    for (size_t j = 0; j < local_sz_samples.size() - t; j++) {
      overlap_sum += SpinConfigurationOverlap(local_sz_samples[j], local_sz_samples[j + t]);
    }
    res[t] = (double) overlap_sum / (local_sz_samples.size() - t) / N - 0.25;
  }
  return res;
}

void PrintProgressBar(int progress, int total) {
  int bar_width = 70; // width of the progress bar

  std::cout << "[";
  int pos = bar_width * progress / total;
  for (int i = 0; i < bar_width; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0 / total) << " %" << std::endl;
}

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

  void ReplicaTest(void); // for check the ergodicity

  void LoadTenData(void);

  void LoadTenData(const std::string &tps_path);

  void DumpData();

  void DumpData(const std::string &tps_path);


  VMCOptimizePara optimize_para;

 private:
  void ReserveSamplesDataSpace_();

  void PrintExecutorInfo_(void);

  void Measure_(void);

  void MCSweep_(size_t &tri_flip_times, size_t &bond_flip_times);

  void WarmUp_(void);

  void MeasureSample_(void);

  void GatherStatistic_(void);

  void SynchronizeConfiguration_(const size_t root = 0); //for the replica test

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
  std::vector<std::vector<TenElemT>> bond_energy_samples_;
  std::vector<size_t> center_configs_; //used to analyze the auto correlation.
  std::vector<std::vector<bool>> local_sz_samples_; // outside is the sample index, inner side is the lattice index.
  // the lattice site number = Lx * Ly * 3,  first the unit cell, then column idx, then row index.

  struct Result {
    TenElemT energy;
    TenElemT en_err;
    std::vector<TenElemT> bond_energys;
    std::vector<double> sz;
    std::vector<TenElemT> energy_auto_corr;
    std::vector<double> spin_auto_corr;
  };
  Result res;
};//KagomeMeasurementExecutor


template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::ReplicaTest() {
  SynchronizeConfiguration_();
  std::vector<double> overlaps;
  overlaps.reserve(optimize_para.mc_samples);
  std::cout << "Random number from worker " << world_.rank() << " : " << u_double_(random_engine) << std::endl;
  const size_t bond_flip_num = lx_ * (ly_ - 1) + (ly_ - 1) * lx_;
  size_t bond_accept_times, tri_accept_times;
  for (size_t sweep = 0; sweep < optimize_para.mc_samples; sweep++) {
    MCSweep_(tri_accept_times, bond_accept_times);
    double bond_accept_ratio = (double) bond_accept_times / double(bond_flip_num);
    double tri_accept_ratio = (double) tri_accept_times / double(bond_flip_num);
    // send-recv configuration
    Configuration config2(ly_, lx_);
    size_t dest = (world_.rank() + 1) % world_.size();
    size_t source = (world_.rank() + world_.size() - 1) % world_.size();
    MPI_Status status;
    int err_msg = MPI_Sendrecv(tps_sample_.config, dest, dest, config2, source, world_.rank(), MPI_Comm(world_),
                               &status);

    // calculate overlap
    double overlap = SpinConfigurationOverlap2(KagomeConfig2Sz(tps_sample_.config), KagomeConfig2Sz(config2));
    overlaps.push_back(overlap);
    if (world_.rank() == kMasterProc && (sweep + 1) % (optimize_para.mc_samples / 10) == 0) {
      PrintProgressBar((sweep + 1), optimize_para.mc_samples);
      std::cout << "Accept Ratio : " << bond_accept_ratio << " , " << tri_accept_ratio << std::endl;
    }
  }
  //DumpData
  std::string replica_overlap_path = "replica_overlap/";
  if (world_.rank() == kMasterProc)
    if (!IsPathExist(replica_overlap_path)) {
      CreatPath(replica_overlap_path);
    }
  world_.barrier();
  DumpVecData(replica_overlap_path + "/replica_overlap" + std::to_string(world_.rank()), overlaps);
}

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
  std::vector<double> bond_energy;
  energy = measurement_solver_.SampleMeasure(&split_index_tps_, &tps_sample_, local_sz, bond_energy);
  energy_samples_.push_back(energy);
  local_sz_samples_.push_back(local_sz);
  bond_energy_samples_.push_back(bond_energy);
  //add more measurement here and the definition of measurement solver
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::GatherStatistic_() {
  TenElemT en_thread = Mean(energy_samples_);
  std::vector<TenElemT> en_list;
  boost::mpi::gather(world_, en_thread, en_list, kMasterProc);
  if (world_.rank() == 0) {
    res.energy = Mean(en_list);
    res.en_err = StandardError(en_list, res.energy);
  }
  res.energy_auto_corr = CalAutoCorrelation(energy_samples_, en_thread);
  const size_t N = 3 * lx_ * ly_; //site number
  std::vector<size_t> sz_sum_thread(N, 0);
  std::vector<std::vector<size_t>> sz_sum_list(N);
  std::vector<size_t> sz_sum(N, 0);
  for (auto &local_sz: local_sz_samples_) {
    for (size_t i = 0; i < N; i++) {
      sz_sum_thread[i] += local_sz[i];
    }
  }
  for (size_t i = 0; i < N; i++) { // i is site index
    boost::mpi::gather(world_, sz_sum_thread[i], sz_sum_list[i], kMasterProc);
    if (world_.rank() == kMasterProc) {
      for (size_t summation: sz_sum_list[i]) { //for every thread's summation result
        sz_sum[i] += summation;
      }
    }
  }
  if (world_.rank() == kMasterProc) {
    res.sz = std::vector<double>(N, 0.0);
    for (size_t i = 0; i < N; i++) {
      res.sz[i] = (double) sz_sum[i] / (double) (optimize_para.mc_samples * world_.size()) - 0.5;
    }
  }
  res.spin_auto_corr = CalSpinAutoCorrelation(local_sz_samples_);

  const size_t bond_num = bond_energy_samples_[0].size();
  res.bond_energys = std::vector<TenElemT>(bond_num);
  std::vector<TenElemT> bond_energy_thread = AveListOfData(bond_energy_samples_);
  for (size_t bond = 0; bond < bond_num; bond++) {
    std::vector<TenElemT> bond_energy_proc_list; //idx is proc number
    boost::mpi::gather(world_, bond_energy_thread[bond], bond_energy_proc_list, kMasterProc);
    if (world_.rank() == kMasterProc) {
      res.bond_energys[bond] = Mean(bond_energy_proc_list);
    }
  }
  std::cout << "Statistic data finished." << std::endl;
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

  if (world_.rank() == kMasterProc) {
    std::ofstream ofs("statistic_summary", std::ofstream::binary);
    ofs.write((const char *) &res.energy, 1 * sizeof(TenElemT));
    ofs.write((const char *) &res.en_err, 1 * sizeof(TenElemT));
    ofs.write((const char *) res.energy_auto_corr.data(), res.energy_auto_corr.size() * sizeof(TenElemT));
    ofs.write((const char *) res.bond_energys.data(), res.bond_energys.size() * sizeof(TenElemT));
    ofs.write((const char *) res.sz.data(), res.sz.size() * sizeof(double));
    ofs.write((const char *) res.spin_auto_corr.data(), res.spin_auto_corr.size() * sizeof(double));
    ofs << std::endl;
    ofs.close();
  }
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
    measurement_solver_(solver), warm_up_(false) {
  TPSSample<TenElemT, QNT>::trun_para = TruncatePara(optimize_para);
  random_engine.seed(std::random_device{}() + world.rank() * 10086);
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
    size_t tri_accept_times, bond_accept_times;
    for (size_t sweep = 0; sweep < optimize_para.mc_warm_up_sweeps; sweep++) {
      MCSweep_(tri_accept_times, bond_accept_times);
    }
    double elasp_time = warm_up_timer.Elapsed();
    std::cout << "Proc " << std::setw(4) << world_.rank() << " warm-up completes T = " << elasp_time << "s."
              << std::endl;
    warm_up_ = true;
  }
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::SynchronizeConfiguration_(const size_t root) {
  Configuration config(tps_sample_.config);
  MPI_BCast(config, root, MPI_Comm(world_));
  if (world_.rank() != root) {
    tps_sample_ = TPSSample<TenElemT, QNT>(split_index_tps_, config);
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
  size_t bond_accept_total_num = 0, tri_accept_total_num = 0;
  size_t bond_num = lx_ * (ly_ - 1) + ly_ * (lx_ - 1);
  size_t tri_num = lx_ * ly_;
  size_t tri_accept_times, bond_accept_times;
  for (size_t sweep = 0; sweep < optimize_para.mc_samples; sweep++) {
    MCSweep_(tri_accept_times, bond_accept_times);
    bond_accept_total_num += bond_accept_times;
    tri_accept_total_num += tri_accept_times;
    MeasureSample_();
    if (world_.rank() == kMasterProc && (sweep + 1) % (optimize_para.mc_samples / 10) == 0) {
      PrintProgressBar((sweep + 1), optimize_para.mc_samples);
    }
  }
  double bond_accept_rate = double(bond_accept_total_num) / double(bond_num * optimize_para.mc_samples);
  double tri_accept_rate = double(tri_accept_times) / double(2 * tri_num * optimize_para.mc_samples);
  std::cout << "Accept Ratio : " << bond_accept_rate << " , " << tri_accept_rate << std::endl;
  GatherStatistic_();
}

template<typename TenElemT, typename QNT, typename MeasurementSolver>
void KagomeMeasurementExecutor<TenElemT, QNT, MeasurementSolver>::MCSweep_(
    size_t &tri_flip_times,
    size_t &bond_flip_times
) {
  for (size_t i = 0; i < optimize_para.mc_sweeps_between_sample; i++) {
    tps_sample_.MCCompressedKagomeLatticeLocalUpdateSweep(split_index_tps_, u_double_, tri_flip_times, bond_flip_times);
  }
}
}//gqpeps

#endif //HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
