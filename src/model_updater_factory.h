// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
#define HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H

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
inline void LogSamplerChoice(const heisenberg_params::PhysicalParams &phys) {
  if (!phys.MCRestrictU1) {
    std::cout << "[info] MCRestrictU1=false (no-U1 sampler requested)."
                 " Using U1 sampler temporarily; non-U1 variant will be wired in a later step."
              << std::endl;
  }
}

// ---------------- VMC dispatcher (by model) ----------------
template<typename TenElemT, typename QNT, typename MCUpdaterT>
inline void RunVmcByModel(EnhancedVMCUpdateParams &params,
                          qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                          MPI_Comm comm,
                          int rank) {
  const std::string model_type = params.physical_params.ModelType.empty() ? "SquareHeisenberg"
                                                                          : params.physical_params.ModelType;
  const double j2 = params.physical_params.J2;
  const bool is_pbc = (params.physical_params.BoundaryCondition == qlpeps::BoundaryCondition::Periodic);

  if (model_type == "SquareHeisenberg") {
    if (is_pbc) {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // PBC TRG
      Model solver(j2);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    } else if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModel; // Heisenberg J2=0 (OBC)
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{}, MCUpdaterT{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // Heisenberg J2!=0 (OBC)
      Model solver(j2);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (is_pbc) {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // XY on PBC via XXZ PBC
      Model solver(/*jz1=*/0.0, /*jxy1=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    } else if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModel; // XY J2=0 => jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // XY J2!=0 => jz=0, jxy=1, jz2=0, jxy2=j2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (is_pbc) {
      throw std::invalid_argument("TriangleHeisenberg PBC is not supported in current PEPS backend.");
    }
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SpinOneHalfTriHeisenbergSqrPEPS; // J2=0
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{}, MCUpdaterT{});
    } else {
      using Model = qlpeps::SpinOneHalfTriJ1J2HeisenbergSqrPEPS; // J2!=0
      Model solver(j2);
      (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, solver, MCUpdaterT{});
    }
    return;
  }

  // Fallback: default to SquareHeisenberg semantics
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModel;
    (void) qlpeps::VmcOptimize(params.CreateVMCOptimizerParams(rank), sitps, comm, Model{}, MCUpdaterT{});
  }
}

// ---------------- Measurement dispatcher (by model) ----------------
template<typename TenElemT, typename QNT, typename MCUpdaterT>
inline void RunMeasureByModel(const heisenberg_params::PhysicalParams &phys,
                              const qlpeps::MCMeasurementParams &measurement_params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm) {
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;
  const bool is_pbc = (phys.BoundaryCondition == qlpeps::BoundaryCondition::Periodic);

  if (model_type == "SquareHeisenberg") {
    if (is_pbc) {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC;
      Model solver(j2);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    } else if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModel;
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel;
      Model solver(j2);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (is_pbc) {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelPBC; // XY on PBC via XXZ PBC
      Model solver(/*jz1=*/0.0, /*jxy1=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    } else if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModel; // XY: jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModel; // XY with J2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (is_pbc) {
      throw std::invalid_argument("TriangleHeisenberg PBC is not supported in current PEPS backend.");
    }
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SpinOneHalfTriHeisenbergSqrPEPS;
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, Model{});
    } else {
      using Model = qlpeps::SpinOneHalfTriJ1J2HeisenbergSqrPEPS;
      Model solver(j2);
      (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, solver);
    }
    return;
  }

  // Fallback
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModel;
    (void) qlpeps::MonteCarloMeasure<TenElemT, QNT, MCUpdaterT, Model>(sitps, measurement_params, comm, Model{});
  }
}

#endif // HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
