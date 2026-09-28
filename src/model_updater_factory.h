// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
#define HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H

#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include "enhanced_params_parser.h"
#include "qlpeps/state/spin_inversion_metadata.h"
#include "qlpeps/vmc_basic/spin_inversion_projected_sample_obc.h"
#include "qlpeps/vmc_basic/mc_updaters/square_nn_spin_inversion_updater_obc.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_inversion_square_xxz_obc.h"
#include "qlpeps/qlpeps.h"

// We keep the actual call sites in the drivers (e.g., square_vmc_update.cpp),
// and use this helper only for detection/logging to keep dependencies localized.
inline void LogSamplerChoice(const heisenberg_params::MonteCarloNumericalParams &mc) {
  if (!mc.MCRestrictU1) {
    std::cout << "[info] MCRestrictU1=false (no-U1 sampler requested)."
                 " Using U1 sampler temporarily; non-U1 variant will be wired in a later step."
              << std::endl;
  }
}

namespace heisenberg_vmcpeps::detail {

template<typename TenElemT,
         typename QNT,
         typename MCUpdaterT,
         typename EnergySolverT,
         template<typename, typename> class ContractorT = qlpeps::BMPSContractor>
inline void ExecuteVmc_(const qlpeps::VMCPEPSOptimizerParams &opt_params,
                        const qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                        const MPI_Comm &comm,
                        const EnergySolverT &solver) {
  using ExecT = qlpeps::VMCPEPSOptimizer<TenElemT, QNT, MCUpdaterT, EnergySolverT, ContractorT>;
  ExecT executor(opt_params, sitps, comm, solver, MCUpdaterT{});
  executor.Execute();
}

template<typename TenElemT,
         typename QNT,
         typename MCUpdaterT,
         typename MeasurementSolverT,
         template<typename, typename> class ContractorT = qlpeps::BMPSContractor>
inline void ExecuteMeasure_(const qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                            const qlpeps::MCMeasurementParams &measurement_params,
                            const MPI_Comm &comm,
                            const MeasurementSolverT &solver) {
  using MeasT = qlpeps::MCPEPSMeasurer<TenElemT, QNT, MCUpdaterT, MeasurementSolverT, ContractorT>;
  MeasT measurer(sitps, measurement_params, comm, solver, MCUpdaterT{});
  measurer.Execute();
}

template<typename TenElemT, typename QNT>
inline void RunVmcByModelOBC_(EnhancedVMCUpdateParams &params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm,
                              int rank) {
  using MCUpdaterT = qlpeps::MCUpdateSquareTNN3SiteExchangeOBC;
  const std::string model_type = params.physical_params.ModelType.empty() ? "SquareHeisenberg"
                                                                          : params.physical_params.ModelType;
  const double j2 = params.physical_params.J2;

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // Heisenberg J2=0 (OBC)
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // Heisenberg J2!=0 (OBC)
      Model solver(j2);
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // XY J2=0 => jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // XY J2!=0 => jz=0, jxy=1, jz2=0, jxy2=j2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::TriangularSpinOneHalfHeisenbergModelOBC; // J2=0
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{});
    } else {
      using Model = qlpeps::TriangularSpinOneHalfJ1J2HeisenbergModelOBC; // J2!=0
      Model solver(j2);
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    }
    return;
  }

  // Fallback: default to SquareHeisenberg semantics
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
    heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{});
  }
}

/** @brief PBC model dispatch for one fixed contraction backend (TRG or HOTRG). */
template<typename TenElemT,
         typename QNT,
         template<typename, typename> class ContractorT>
inline void RunVmcByModelPBCWith_(const qlpeps::VMCPEPSOptimizerParams &opt_params,
                                  const heisenberg_params::PhysicalParams &phys,
                                  qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                  MPI_Comm comm) {
  using MCUpdaterT = qlpeps::MCUpdateSquareNNExchangePBC;
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;

  if (model_type == "TriangleHeisenberg") {
    throw std::invalid_argument("TriangleHeisenberg PBC is not supported in current PEPS backend.");
  }

  if (model_type == "SquareXY") {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // XY on PBC via XXZ PBC
    Model solver(/*jz1=*/0.0, /*jxy1=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
    heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT, Model, ContractorT>(
        opt_params, sitps, comm, solver);
    return;
  }

  // Default: SquareHeisenberg semantics on PBC.
  {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC;
    Model solver(j2);
    heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT, Model, ContractorT>(
        opt_params, sitps, comm, solver);
  }
}

/** @brief Run PBC VMC, choosing the backend held by the parsed ContractorParams. */
template<typename TenElemT, typename QNT>
inline void RunVmcByModelPBC_(EnhancedVMCUpdateParams &params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm,
                              int rank) {
  const qlpeps::VMCPEPSOptimizerParams opt_params = params.CreateVMCOptimizerParams(rank);
  if (opt_params.contractor_params.IsHOTRG()) {
    RunVmcByModelPBCWith_<TenElemT, QNT, qlpeps::HOTRGContractor>(
        opt_params, params.physical_params, sitps, comm);
    return;
  }
  RunVmcByModelPBCWith_<TenElemT, QNT, qlpeps::TRGContractor>(
      opt_params, params.physical_params, sitps, comm);
}

/** @brief Optimize the coherent spin-inversion state with shared PEPS parameters. */
template<typename TenElemT, typename QNT, int Parity>
inline void RunSpinInversionVmc_(EnhancedVMCUpdateParams &params,
                                 const qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                 MPI_Comm comm, int rank) {
  using Sample = qlpeps::SpinInversionProjectedSampleOBC<TenElemT, QNT, Parity>;
  using Updater = qlpeps::MCUpdateSquareNNSpinInversionExchangeOBC<Parity>;
  using Model = qlpeps::SpinInversionSquareXXZModelOBC;
  using Executor = qlpeps::VMCPEPSOptimizer<TenElemT, QNT, Updater, Model,
                                          qlpeps::BMPSContractor, Sample>;
  auto optimizer_params = params.CreateVMCOptimizerParams(rank);
  const auto &config = optimizer_params.mc_params.initial_config;
  size_t up = 0, down = 0;
  for (size_t row = 0; row < config.rows(); ++row) {
    for (size_t col = 0; col < config.cols(); ++col) {
      up += config({row, col}) == 0;
      down += config({row, col}) == 1;
    }
  }
  int invalid = (up != down || up + down != config.rows() * config.cols());
  int any_invalid = 0;
  ::MPI_Allreduce(&invalid, &any_invalid, 1, MPI_INT, MPI_MAX, comm);
  if (any_invalid) throw std::invalid_argument("Spin inversion requires Sz=0 configurations on every rank.");
  optimizer_params.mc_params.assume_initial_config_thermalized = false;
  optimizer_params.tps_dump_base_name = params.io_params.wavefunction_base;
  Model model(params.physical_params.ModelType == "SquareXY" ? 0.0 : 1.0, 1.0);
  Executor executor(optimizer_params, sitps, comm, model, Updater{});
  if (rank == 0) {
    std::cout << "SpinInversionParity=" << Parity
              << ": optimizing psi(x) + parity*psi(Fx); configured warm-up is always run.\n"
              << "Saved tensors require this parity; plain mc_measure is unsupported."
              << std::endl;
  }
  executor.Execute();
}

// ---------------- Measurement dispatcher (by model) ----------------
template<typename TenElemT, typename QNT>
inline void RunMeasureByModelOBC_(const heisenberg_params::PhysicalParams &phys,
                                  const qlpeps::MCMeasurementParams &measurement_params,
                                  qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                  MPI_Comm comm) {
  using MCUpdaterT = qlpeps::MCUpdateSquareTNN3SiteExchangeOBC;
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC;
      Model solver(j2);
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // XY: jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, solver);
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // XY with J2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::TriangularSpinOneHalfHeisenbergModelOBC;
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::TriangularSpinOneHalfJ1J2HeisenbergModelOBC;
      Model solver(j2);
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  // Fallback
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
    heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, Model{});
  }
}

/** @brief PBC measurement dispatch for one fixed contraction backend (TRG or HOTRG). */
template<typename TenElemT,
         typename QNT,
         template<typename, typename> class ContractorT>
inline void RunMeasureByModelPBCWith_(const heisenberg_params::PhysicalParams &phys,
                                      const qlpeps::MCMeasurementParams &measurement_params,
                                      qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                      MPI_Comm comm) {
  using MCUpdaterT = qlpeps::MCUpdateSquareNNExchangePBC;
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;

  if (model_type == "TriangleHeisenberg") {
    throw std::invalid_argument("TriangleHeisenberg PBC is not supported in current PEPS backend.");
  }

  if (model_type == "SquareXY") {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // XY on PBC via XXZ PBC
    Model solver(/*jz1=*/0.0, /*jxy1=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
    heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT, Model, ContractorT>(
        sitps, measurement_params, comm, solver);
    return;
  }

  // Default: SquareHeisenberg semantics on PBC.
  {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC;
    Model solver(j2);
    heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT, Model, ContractorT>(
        sitps, measurement_params, comm, solver);
  }
}

/** @brief Run PBC measurement, choosing the backend held by the parsed ContractorParams. */
template<typename TenElemT, typename QNT>
inline void RunMeasureByModelPBC_(const heisenberg_params::PhysicalParams &phys,
                                  const qlpeps::MCMeasurementParams &measurement_params,
                                  qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                  MPI_Comm comm) {
  if (measurement_params.contractor_params.IsHOTRG()) {
    RunMeasureByModelPBCWith_<TenElemT, QNT, qlpeps::HOTRGContractor>(
        phys, measurement_params, sitps, comm);
    return;
  }
  RunMeasureByModelPBCWith_<TenElemT, QNT, qlpeps::TRGContractor>(
      phys, measurement_params, sitps, comm);
}

} // namespace heisenberg_vmcpeps::detail

// ---------------- Public dispatchers ----------------
template<typename TenElemT, typename QNT>
inline void RunVmcByModel(EnhancedVMCUpdateParams &params,
                          qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                          MPI_Comm comm,
                          int rank) {
  if (params.spin_inversion_parity == 1) {
    return heisenberg_vmcpeps::detail::RunSpinInversionVmc_<TenElemT, QNT, 1>(params, sitps, comm, rank);
  }
  if (params.spin_inversion_parity == -1) {
    return heisenberg_vmcpeps::detail::RunSpinInversionVmc_<TenElemT, QNT, -1>(params, sitps, comm, rank);
  }
  const bool is_pbc = (params.physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic);
  if (is_pbc) {
    heisenberg_vmcpeps::detail::RunVmcByModelPBC_<TenElemT, QNT>(params, sitps, comm, rank);
  } else {
    heisenberg_vmcpeps::detail::RunVmcByModelOBC_<TenElemT, QNT>(params, sitps, comm, rank);
  }
}

template<typename TenElemT, typename QNT>
inline void RunMeasureByModel(const heisenberg_params::PhysicalParams &phys,
                              const qlpeps::MCMeasurementParams &measurement_params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm) {
  const bool is_pbc = (phys.BoundaryCondition == qlpeps::BoundaryCondition::Periodic);
  if (is_pbc) {
    heisenberg_vmcpeps::detail::RunMeasureByModelPBC_<TenElemT, QNT>(phys, measurement_params, sitps, comm);
  } else {
    heisenberg_vmcpeps::detail::RunMeasureByModelOBC_<TenElemT, QNT>(phys, measurement_params, sitps, comm);
  }
}

#endif // HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
