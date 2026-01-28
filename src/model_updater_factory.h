// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
#define HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H

#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include "enhanced_params_parser.h"
#include "qlpeps/qlpeps.h"
// Future: include no-U1 updater types when available

enum class ModelKind {
  SquareHeisenberg,
  SquareXY,
  TriangleHeisenberg,
  Unknown
};

inline ModelKind DetectModelKind(const std::string &model_type) {
  if (model_type == "SquareHeisenberg") return ModelKind::SquareHeisenberg;
  if (model_type == "SquareXY") return ModelKind::SquareXY;
  if (model_type == "TriangleHeisenberg") return ModelKind::TriangleHeisenberg;
  return ModelKind::Unknown;
}

// Current minimal mapping to avoid code bloat.
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
  using MCUpdaterT = qlpeps::MCUpdateSquareTNN3SiteExchange;
  const std::string model_type = params.physical_params.ModelType.empty() ? "SquareHeisenberg"
                                                                          : params.physical_params.ModelType;
  const double j2 = params.physical_params.J2;

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // Heisenberg J2=0 (OBC)
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // Heisenberg J2!=0 (OBC)
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
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // XY J2!=0 => jz=0, jxy=1, jz2=0, jxy2=j2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SpinOneHalfTriHeisenbergSqrPEPS; // J2=0
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT>(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{});
    } else {
      using Model = qlpeps::SpinOneHalfTriJ1J2HeisenbergSqrPEPS; // J2!=0
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

template<typename TenElemT, typename QNT>
inline void RunVmcByModelPBC_(EnhancedVMCUpdateParams &params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm,
                              int rank) {
  using MCUpdaterT = qlpeps::MCUpdateSquareNNExchangePBC;
  const std::string model_type = params.physical_params.ModelType.empty() ? "SquareHeisenberg"
                                                                          : params.physical_params.ModelType;
  const double j2 = params.physical_params.J2;

  if (model_type == "TriangleHeisenberg") {
    throw std::invalid_argument("TriangleHeisenberg PBC is not supported in current PEPS backend.");
  }

  if (model_type == "SquareXY") {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // XY on PBC via XXZ PBC
    Model solver(/*jz1=*/0.0, /*jxy1=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
    heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT, Model, qlpeps::TRGContractor>(
        params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
    return;
  }

  // Default: SquareHeisenberg semantics on PBC/TRG.
  {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // PBC TRG
    Model solver(j2);
    heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, MCUpdaterT, Model, qlpeps::TRGContractor>(
        params.CreateVMCOptimizerParams(rank), sitps, comm, solver);
  }
}

// ---------------- Measurement dispatcher (by model) ----------------
template<typename TenElemT, typename QNT>
inline void RunMeasureByModelOBC_(const heisenberg_params::PhysicalParams &phys,
                                  const qlpeps::MCMeasurementParams &measurement_params,
                                  qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                  MPI_Comm comm) {
  using MCUpdaterT = qlpeps::MCUpdateSquareTNN3SiteExchange;
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel;
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
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // XY with J2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SpinOneHalfTriHeisenbergSqrPEPS;
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::SpinOneHalfTriJ1J2HeisenbergSqrPEPS;
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

template<typename TenElemT, typename QNT>
inline void RunMeasureByModelPBC_(const heisenberg_params::PhysicalParams &phys,
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
    heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT, Model, qlpeps::TRGContractor>(
        sitps, measurement_params, comm, solver);
    return;
  }

  // Default: SquareHeisenberg semantics on PBC/TRG.
  {
    using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC;
    Model solver(j2);
    heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, MCUpdaterT, Model, qlpeps::TRGContractor>(
        sitps, measurement_params, comm, solver);
  }
}

} // namespace heisenberg_vmcpeps::detail

// ---------------- Public dispatchers ----------------
template<typename TenElemT, typename QNT>
inline void RunVmcByModel(EnhancedVMCUpdateParams &params,
                          qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                          MPI_Comm comm,
                          int rank) {
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
