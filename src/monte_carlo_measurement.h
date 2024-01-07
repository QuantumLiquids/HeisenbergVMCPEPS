//
// Created by haoxinwang on 02/11/2023.
//

#ifndef HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
#define HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H

#include "gqpeps/algorithm/vmc_update/vmc_peps.h"

namespace gqpeps {
using namespace gqten;

std::vector<bool> KagomeConfig2Sz(
    const Configuration &config
) {
  std::vector<bool> local_sz;
  local_sz.reserve(config.size() * 3 - config.rows() - config.cols());
  for (size_t row = 0; row < config.rows(); row++) {
    for (size_t col = 0; col < config.cols(); col++) {
      size_t local_config = config({row, col});
      local_sz.push_back(local_config & 1); //left upper site
      if (row < config.cols() - 1)    // remove corner
        local_sz.push_back(local_config >> 1 & 1);//lower site
      if (col < config.cols() - 1)  //remove corner
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

/**
 * Calculate the spin auto-correlation from the local_sz_samples
 * The auto-correlation is defined as
 *
 * 1/N Sum_i <S_i(t)*S_i(t+delta t)> - 1/N Sum_i <S_i>^2
 *
 * Where S_i = 0 or 1; or equivalently \pm 0.5
 * and correspondingly we assume <S_i> = 0.5; or 0. But this assumption may not correct.
 * @param local_sz_samples
 * @return
 */
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

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
class MonteCarloMeasurementExecutor : public Executor {
 public:
  using Tensor = GQTensor<TenElemT, QNT>;
  using TPST = TPS<TenElemT, QNT>;
  using SITPST = SplitIndexTPS<TenElemT, QNT>;
  using IndexT = Index<QNT>;

  //Load Data from path
  MonteCarloMeasurementExecutor(const VMCOptimizePara &optimize_para,
                                const size_t ly, const size_t lx,
                                const boost::mpi::communicator &world,
                                const MeasurementSolver &solver = MeasurementSolver());

  MonteCarloMeasurementExecutor(const VMCOptimizePara &optimize_para,
                                const SITPST &sitpst_init,
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

  std::vector<double> MCSweep_(void);

  void WarmUp_(void);

  void MeasureSample_(void);

  void GatherStatistic_(void);

  void SynchronizeConfiguration_(const size_t root = 0); //for the replica test

  boost::mpi::communicator world_;

  size_t lx_; //cols
  size_t ly_; //rows

  SITPST split_index_tps_;

  WaveFunctionComponentType tps_sample_;

  std::uniform_real_distribution<double> u_double_;

  bool warm_up_;

  MeasurementSolver measurement_solver_;

  // observable
  struct SampleData {
    std::vector<TenElemT> energy_samples;
    std::vector<std::vector<TenElemT>> bond_energy_samples;
    std::vector<std::vector<TenElemT>>
        one_point_function_samples; // outside is the sample index, inner side is the lattice index.
    std::vector<std::vector<TenElemT>>
        two_point_function_samples;

    void Reserve(const size_t sample_num) {
      energy_samples.reserve(sample_num);
      bond_energy_samples.reserve(sample_num);
      one_point_function_samples.reserve(sample_num);
      two_point_function_samples.reserve(sample_num);
    }
  } sample_data_;
  // the lattice site number = Lx * Ly * 3,  first the unit cell, then column idx, then row index.

  struct Result {
    TenElemT energy;
    TenElemT en_err;
    std::vector<TenElemT> energy_auto_corr;
    std::vector<TenElemT> bond_energys;
    std::vector<double> one_point_functions;
    std::vector<double> two_point_functions;
    std::vector<double> one_auto_corr;
  };
  Result res;
};//MonteCarloMeasurementExecutor


template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::ReplicaTest() {
  SynchronizeConfiguration_();
  std::vector<double> overlaps;
  overlaps.reserve(optimize_para.mc_samples);
//  std::cout << "Random number from worker " << world_.rank() << " : " << u_double_(random_engine) << std::endl;
  std::vector<double> accept_rates_accum;
  for (size_t sweep = 0; sweep < optimize_para.mc_samples; sweep++) {
    std::vector<double> accept_rates = MCSweep_();
    if (sweep == 0) {
      accept_rates_accum = accept_rates;
    } else {
      for (size_t i = 0; i < accept_rates_accum.size(); i++) {
        accept_rates_accum[i] += accept_rates[i];
      }
    }
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

      auto accept_rates_avg = accept_rates_accum;
      for (double &rates : accept_rates_avg) {
        rates /= double(sweep + 1);
      }
      std::cout << "Accept rate = [";
      for (double &rate : accept_rates_avg) {
        std::cout << std::setw(5) << std::fixed << std::setprecision(2) << rate;
      }
      std::cout << "]";
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
  tps_sample_.config.Dump(optimize_para.wavefunction_path, world_.rank());
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT,
                                   QNT,
                                   WaveFunctionComponentType,
                                   MeasurementSolver>::ReserveSamplesDataSpace_(
    void) {
  sample_data_.Reserve(optimize_para.mc_samples);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::MeasureSample_() {
  TenElemT energy;
  std::vector<bool> local_sz;
  std::vector<double> bond_energy;
  energy = measurement_solver_.SampleMeasure(&split_index_tps_, &tps_sample_, local_sz, bond_energy);
  energy_samples.push_back(energy);
  bond_energy_samples.push_back(bond_energy);
  one_point_function_samples.push_back(local_sz);
  //add more measurement here and the definition of measurement solver
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::GatherStatistic_() {
  TenElemT en_thread = Mean(energy_samples_);
  std::vector<TenElemT> en_list;
  boost::mpi::gather(world_, en_thread, en_list, kMasterProc);
  if (world_.rank() == 0) {
    res.energy = Mean(en_list);
    res.en_err = StandardError(en_list, res.energy);
  }
  res.energy_auto_corr = CalAutoCorrelation(energy_samples_, en_thread);
  const size_t N = local_sz_samples_[0].size();
  //physical site number.
  //3 * lx * ly for complete unit cell lattice; 3 * lx * ly - lx - ly for smooth boundary lattice.
  std::vector<size_t> sz_sum_thread(N, 0);
  std::vector<std::vector<size_t>> sz_sum_list(N);
  std::vector<size_t> sz_sum(N, 0);
  for (auto &local_sz : local_sz_samples_) {
    for (size_t i = 0; i < N; i++) {
      sz_sum_thread[i] += local_sz[i];
    }
  }
  for (size_t i = 0; i < N; i++) { // i is site index
    boost::mpi::gather(world_, sz_sum_thread[i], sz_sum_list[i], kMasterProc);
    if (world_.rank() == kMasterProc) {
      for (size_t summation : sz_sum_list[i]) { //for every thread's summation result
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
  std::cout << "Rank " << world_.rank() << ": statistic data finished." << std::endl;
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void
MonteCarloMeasurementExecutor<TenElemT,
                              QNT,
                              WaveFunctionComponentType,
                              MeasurementSolver>::DumpData(const std::string &tps_path) {
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

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
MonteCarloMeasurementExecutor<TenElemT,
                              QNT,
                              WaveFunctionComponentType,
                              MeasurementSolver>::MonteCarloMeasurementExecutor(
    const VMCOptimizePara &optimize_para,
    const size_t ly, const size_t lx,
    const boost::mpi::communicator &world,
    const MeasurementSolver &solver):
    optimize_para(optimize_para), world_(world), lx_(lx), ly_(ly),
    split_index_tps_(ly, lx), tps_sample_(ly, lx),
    u_double_(0, 1), warm_up_(false),
    measurement_solver_(solver) {
  WaveFunctionComponentType::trun_para = BMPSTruncatePara(optimize_para);
  random_engine.seed(std::random_device{}() + world.rank() * 10086);
  LoadTenData();
  ReserveSamplesDataSpace_();
  PrintExecutorInfo_();
  this->SetStatus(ExecutorStatus::INITED);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
MonteCarloMeasurementExecutor<TenElemT,
                              QNT,
                              WaveFunctionComponentType,
                              MeasurementSolver>::MonteCarloMeasurementExecutor(
    const VMCOptimizePara &optimize_para,
    const SITPST &sitpst_init,
    const boost::mpi::communicator &world,
    const MeasurementSolver &solver):
    optimize_para(optimize_para), world_(world),
    lx_(sitpst_init.cols()),
    ly_(sitpst_init.rows()),
    split_index_tps_(sitpst_init),
    tps_sample_(ly_, lx_),
    u_double_(0, 1), warm_up_(false),
    measurement_solver_(solver) {
  WaveFunctionComponentType::trun_para = BMPSTruncatePara(optimize_para);
  random_engine.seed(std::random_device{}() + world.rank() * 10086);
  tps_sample_ = WaveFunctionComponentType(sitpst_init, optimize_para.init_config);
  ReserveSamplesDataSpace_();
  PrintExecutorInfo_();
  this->SetStatus(ExecutorStatus::INITED);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT,
                                   QNT,
                                   WaveFunctionComponentType,
                                   MeasurementSolver>::PrintExecutorInfo_(void) {
  if (world_.rank() == kMasterProc) {
    std::cout << std::left;  // Set left alignment for the output
    std::cout << "\n";
    std::cout << "=====> MONTE-CARLO MEASUREMENT PROGRAM FOR PEPS <=====" << "\n";
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

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::Execute(void) {
  SetStatus(ExecutorStatus::EXEING);
  Measure_();
  DumpData();
  SetStatus(ExecutorStatus::FINISH);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::WarmUp_(void) {
  if (!warm_up_) {
    Timer warm_up_timer("warm_up");
    for (size_t sweep = 0; sweep < optimize_para.mc_warm_up_sweeps; sweep++) {
      auto accept_rates = MCSweep_();
    }
    double elasp_time = warm_up_timer.Elapsed();
    std::cout << "Proc " << std::setw(4) << world_.rank() << " warm-up completes T = " << elasp_time << "s."
              << std::endl;
    warm_up_ = true;
  }
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT,
                                   QNT,
                                   WaveFunctionComponentType,
                                   MeasurementSolver>::SynchronizeConfiguration_(
    const size_t root) {
  Configuration config(tps_sample_.config);
  MPI_BCast(config, root, MPI_Comm(world_));
  if (world_.rank() != root) {
    tps_sample_ = WaveFunctionComponentType(split_index_tps_, config);
  }
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::LoadTenData(void) {
  LoadTenData(optimize_para.wavefunction_path);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT,
                                   QNT,
                                   WaveFunctionComponentType,
                                   MeasurementSolver>::LoadTenData(const std::string &tps_path) {
  if (!split_index_tps_.Load(tps_path)) {
    std::cout << "Loading TPS files fails." << std::endl;
    exit(-1);
  }
  Configuration config(ly_, lx_);
  bool load_config = config.Load(tps_path, world_.rank());
  if (load_config) {
    tps_sample_ = WaveFunctionComponentType(split_index_tps_, config);
  } else {
    std::cout << "Loading configuration in rank " << world_.rank()
              << " fails. Random generate it and warm up."
              << std::endl;
    tps_sample_ = WaveFunctionComponentType(split_index_tps_, optimize_para.init_config);
    WarmUp_();
  }
  warm_up_ = true;
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::DumpData(void) {
  DumpData(optimize_para.wavefunction_path);
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
void MonteCarloMeasurementExecutor<TenElemT, QNT, WaveFunctionComponentType, MeasurementSolver>::Measure_(void) {
  std::vector<double> accept_rates_accum;
  for (size_t sweep = 0; sweep < optimize_para.mc_samples; sweep++) {
    std::vector<double> accept_rates = MCSweep_();
    if (sweep == 0) {
      accept_rates_accum = accept_rates;
    } else {
      for (size_t i = 0; i < accept_rates_accum.size(); i++) {
        accept_rates_accum[i] += accept_rates[i];
      }
    }
    MeasureSample_();
    if (world_.rank() == kMasterProc && (sweep + 1) % (optimize_para.mc_samples / 10) == 0) {
      PrintProgressBar((sweep + 1), optimize_para.mc_samples);
    }
  }
  std::vector<double> accept_rates_avg = accept_rates_accum;
  for (double &rates : accept_rates_avg) {
    rates /= double(optimize_para.mc_samples);
  }
  std::cout << "Accept rate = [";
  for (double &rate : accept_rates_avg) {
    std::cout << std::setw(5) << std::fixed << std::setprecision(2) << rate;
  }
  std::cout << "]";
  GatherStatistic_();
}

template<typename TenElemT, typename QNT, typename WaveFunctionComponentType, typename MeasurementSolver>
std::vector<double> MonteCarloMeasurementExecutor<TenElemT,
                                                  QNT,
                                                  WaveFunctionComponentType,
                                                  MeasurementSolver>::MCSweep_(
    void) {
  std::vector<double> accept_rates;
  for (size_t i = 0; i < optimize_para.mc_sweeps_between_sample; i++) {
    tps_sample_.MonteCarloSweepUpdate(split_index_tps_, unit_even_distribution, accept_rates);
  }
  return accept_rates;
}
}//gqpeps

#endif //HEISENBERGVMCPEPS_KAGOME_MC_MEASURE_H
