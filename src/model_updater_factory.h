// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
#define HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H

#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

// The opt-in axis update (VisitOBCUpdater) needs PEPS main at or after 68452fe. PEPS reports
// version 0.2.2 both before and after these headers arrived, so find_package(PEPS 0.2.2) in
// CMakeLists.txt cannot reject an older install; this check names the requirement instead.
#if !__has_include("qlpeps/vmc_basic/mc_updaters/square_axis_autoregressive_updater_obc.h") || \
    !__has_include("qlpeps/vmc_basic/mc_updaters/mc_update_sequence.h")
#error "This PEPS install lacks the row/column axis update (MCUpdateSquareAxisAutoregressiveOBC, MCUpdateSequence). Build against PEPS main at or after 68452fe; see skills/changelog.md."
#endif

#include "enhanced_params_parser.h"
#include "qlpeps/vmc_basic/spin_inversion_metadata.h"
#include "qlpeps/vmc_basic/spin_inversion_projected_sample_obc.h"
#include "qlpeps/vmc_basic/point_group_projected_sample_obc.h"
#include "qlpeps/vmc_basic/mc_updaters/square_nn_point_group_updater_obc.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/point_group_square_xxz_obc.h"
#include "qlpeps/vmc_basic/mc_updaters/square_nn_spin_inversion_updater_obc.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_inversion_square_xxz_obc.h"
#include "qlpeps/qlpeps.h"

/**
 * @brief Map the `MCAxisUpdate*` keys to the PEPS axis-update parameters.
 *
 * - Slice MPS: `BMPSSliceSamplerParams(MCAxisUpdateDmin, MCAxisUpdateDmax, MCAxisUpdateTruncErr)`.
 * - Count table: `LabelCountTable({{1}, {0}})` (N_up; label 0 = up, see qldouble.h) only for dense
 *   tensors (TrivialRepQN) with `MCRestrictU1=true`, so that every axis move keeps the global S_z.
 *   U1QN tensors carry S_z themselves, so a table would be redundant and none is set. With
 *   `MCRestrictU1=false` and dense tensors none is set either: the axis moves then change S_z
 *   while the local 3-site exchange conserves it, so the combined chain samples all S_z sectors.
 * - Axes: rows and columns (the library default); no exactness self-check.
 *
 * @pre @p axis was validated (AxisUpdateParams::Validate); the library throws
 *      std::invalid_argument for slice-MPS parameters out of range.
 */
template<typename QNT>
inline qlpeps::AxisAutoregressiveUpdateParams MakeAxisAutoregressiveUpdateParams(
    const heisenberg_params::AxisUpdateParams &axis, bool mc_restrict_u1) {
  qlpeps::AxisAutoregressiveUpdateParams params(qlpeps::BMPSSliceSamplerParams(
      axis.MCAxisUpdateDmin, axis.MCAxisUpdateDmax, axis.MCAxisUpdateTruncErr));
  if (mc_restrict_u1 && std::is_same_v<QNT, qlten::special_qn::TrivialRepQN>) {
    params.count_constraint = qlpeps::LabelCountTable({{1}, {0}});
  }
  return params;
}

/**
 * @brief Call @p f with the OBC sweep updater the parameters select.
 *
 * - `MCAxisUpdate` false (default): `f(qlpeps::MCUpdateSquareTNN3SiteExchangeOBC{})`, the
 *   historical updater, default-constructed as before.
 * - `MCAxisUpdate` true: `f(MCUpdateSequence<MCUpdateSquareAxisAutoregressiveOBC,
 *   MCUpdateSquareTNN3SiteExchangeOBC>)`: per sweep one rejection-free row-and-column pass, then
 *   one local 3-site exchange sweep, with the parameters of MakeAxisAutoregressiveUpdateParams().
 *
 * The two updaters have different types, so @p f must be generic (e.g. `[&](auto updater)`).
 */
template<typename QNT, typename F>
inline void VisitOBCUpdater(const heisenberg_params::AxisUpdateParams &axis,
                            const heisenberg_params::MonteCarloNumericalParams &mc,
                            F &&f) {
  using LocalUpdater = qlpeps::MCUpdateSquareTNN3SiteExchangeOBC;
  if (!axis.MCAxisUpdate) {
    std::forward<F>(f)(LocalUpdater{});
    return;
  }
  using AxisUpdater = qlpeps::MCUpdateSquareAxisAutoregressiveOBC;
  std::forward<F>(f)(qlpeps::MCUpdateSequence<AxisUpdater, LocalUpdater>(
      AxisUpdater(MakeAxisAutoregressiveUpdateParams<QNT>(axis, mc.MCRestrictU1)), LocalUpdater{}));
}

/**
 * @brief Log the sampler choice of the OBC drivers.
 *
 * With `MCAxisUpdate` off the output is unchanged: every rank prints the MCRestrictU1=false
 * notice (the local updater always conserves S_z), nothing otherwise. With `MCAxisUpdate` on,
 * rank 0 prints the composed updater, the slice-MPS parameters, the count-table choice and the
 * meaning of the `[MC acceptance]` components.
 */
template<typename QNT>
inline void LogSamplerChoice(const heisenberg_params::MonteCarloNumericalParams &mc,
                             const heisenberg_params::AxisUpdateParams &axis,
                             int rank) {
  if (!axis.MCAxisUpdate) {
    if (!mc.MCRestrictU1) {
      std::cout << "[info] MCRestrictU1=false (no-U1 sampler requested)."
                   " Using U1 sampler temporarily; non-U1 variant will be wired in a later step."
                << std::endl;
    }
    return;
  }
  if (rank != 0) return;
  const auto params = MakeAxisAutoregressiveUpdateParams<QNT>(axis, mc.MCRestrictU1);
  const char *table = params.count_constraint.has_value()
      ? "N_up {{1},{0}} (dense tensors, MCRestrictU1=true: S_z is conserved)"
      : (std::is_same_v<QNT, qlten::special_qn::TrivialRepQN>
             ? "none (dense tensors, MCRestrictU1=false: axis moves change S_z, all S_z sectors are sampled)"
             : "none (U1QN tensors conserve S_z)");
  std::cout << "[info] MCAxisUpdate=true: each OBC sweep runs MCUpdateSquareAxisAutoregressiveOBC "
               "(rejection-free row and column moves), then MCUpdateSquareTNN3SiteExchangeOBC.\n"
            << "[info]   slice MPS: D_min=" << params.slice_params.D_min
            << " D_max=" << params.slice_params.D_max
            << " trunc_err=" << params.slice_params.trunc_err << "\n"
            << "[info]   count table: " << table << "\n"
            << "[info]   [MC acceptance] components: 0 = changed rows / Ly, 1 = changed columns / Lx, "
               "2 = local 3-site exchange acceptance; [MC updater] lines carry child0.axis.* counts."
            << std::endl;
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
                        const EnergySolverT &solver,
                        MCUpdaterT updater = MCUpdaterT{}) {
  using ExecT = qlpeps::VMCPEPSOptimizer<TenElemT, QNT, MCUpdaterT, EnergySolverT, ContractorT>;
  ExecT executor(opt_params, sitps, comm, solver, std::move(updater));
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
                            const MeasurementSolverT &solver,
                            MCUpdaterT updater = MCUpdaterT{}) {
  using MeasT = qlpeps::MCPEPSMeasurer<TenElemT, QNT, MCUpdaterT, MeasurementSolverT, ContractorT>;
  MeasT measurer(sitps, measurement_params, comm, solver, std::move(updater));
  measurer.Execute();
}

template<typename TenElemT, typename QNT>
inline void RunVmcByModelOBC_(EnhancedVMCUpdateParams &params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm,
                              int rank) {
  const std::string model_type = params.physical_params.ModelType.empty() ? "SquareHeisenberg"
                                                                          : params.physical_params.ModelType;
  const double j2 = params.physical_params.J2;
  // Every model runs the updater VisitOBCUpdater selects: the 3-site exchange alone (default)
  // or the axis update followed by it (MCAxisUpdate=true). The initial configuration is built
  // before the updater, as before.
  const auto run = [&](const auto &solver) {
    const qlpeps::VMCPEPSOptimizerParams opt_params = params.CreateVMCOptimizerParams(rank);
    VisitOBCUpdater<QNT>(params.axis_update_params, params.mc_params, [&](auto updater) {
      heisenberg_vmcpeps::detail::ExecuteVmc_<TenElemT, QNT, decltype(updater)>(
          opt_params, sitps, comm, solver, std::move(updater));
    });
  };

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // Heisenberg J2=0 (OBC)
      run(Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // Heisenberg J2!=0 (OBC)
      Model solver(j2);
      run(solver);
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // XY J2=0 => jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      run(solver);
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // XY J2!=0 => jz=0, jxy=1, jz2=0, jxy2=j2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      run(solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::TriangularSpinOneHalfHeisenbergModelOBC; // J2=0
      run(Model{});
    } else {
      using Model = qlpeps::TriangularSpinOneHalfJ1J2HeisenbergModelOBC; // J2!=0
      Model solver(j2);
      run(solver);
    }
    return;
  }

  // Fallback: default to SquareHeisenberg semantics
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
    run(Model{});
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
              << "Measure saved tensors with the same SpinInversionParity."
              << std::endl;
  }
  executor.Execute();
}

/** @brief Validate every rank's loaded configuration before projected sampling starts. */
inline void RequireProjectedConfiguration_(const qlpeps::Configuration &config,
                                          const qlpeps::PointGroupProjectionParams &projection,
                                          MPI_Comm comm) {
  size_t up = 0, down = 0;
  for (size_t row = 0; row < config.rows(); ++row) {
    for (size_t col = 0; col < config.cols(); ++col) {
      up += config({row, col}) == 0;
      down += config({row, col}) == 1;
    }
  }
  const bool valid = up + down == config.rows() * config.cols() &&
      (projection.spin_inversion_parity == 0 || up == down);
  qlpeps::RequireCollectivelyValid(valid, comm,
      "Projected sampling requires binary configurations and Sz=0 when spin inversion is enabled.");
}

/** @brief XXZ couplings matching the unprojected Heisenberg/XY J1-J2 models. */
inline qlpeps::PointGroupSquareXXZModelOBC MakePointGroupModel_(
    const heisenberg_params::PhysicalParams &physics) {
  const bool xy = physics.ModelType == "SquareXY";
  return qlpeps::PointGroupSquareXXZModelOBC(
      xy ? 0.0 : 1.0, 1.0, xy ? 0.0 : physics.J2, physics.J2);
}

/** @brief Report the coherent ansatz and its current contraction cost. */
inline void LogPointGroupProjection_(const qlpeps::PointGroupProjectionParams &projection,
                                    MPI_Comm comm) {
  int rank = 0;
  ::MPI_Comm_rank(comm, &rank);
  if (rank != 0) return;
  std::cout << "PointGroup=" << projection.group << " PointGroupIrrep=" << projection.irrep
            << " SpinInversionParity=" << projection.spin_inversion_parity
            << ": coherent projected sampling; configured warm-up is always run.\n"
            << "Each proposed configuration recomputes all nonzero symmetry branches; "
               "this reference path is more expensive than the plain cached sampler."
            << std::endl;
  if (projection.group == "D4" && projection.irrep == "E") {
    std::cout << "D4 E projects onto the complete E isotypic subspace; it does not select "
                 "a rotation eigenvector within the doublet." << std::endl;
  }
}

/** @brief Optimize the spatially projected PEPS with one shared set of tensors. */
template<typename TenElemT, typename QNT>
inline void RunPointGroupVmc_(EnhancedVMCUpdateParams &params,
                             const qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                             MPI_Comm comm, int rank) {
  using Sample = qlpeps::PointGroupProjectedSampleOBC<TenElemT, QNT>;
  using Updater = qlpeps::MCUpdateSquareNNPointGroupExchangeOBC;
  using Model = qlpeps::PointGroupSquareXXZModelOBC;
  using Executor = qlpeps::VMCPEPSOptimizer<TenElemT, QNT, Updater, Model,
                                         qlpeps::BMPSContractor, Sample>;
  auto optimizer_params = params.CreateVMCOptimizerParams(rank);
  RequireProjectedConfiguration_(optimizer_params.mc_params.initial_config,
                                 params.point_group_projection, comm);
  optimizer_params.tps_dump_base_name = params.io_params.wavefunction_base;
  const Model model = MakePointGroupModel_(params.physical_params);
  LogPointGroupProjection_(params.point_group_projection, comm);
  Executor executor(optimizer_params, sitps, comm, model, Updater{});
  executor.Execute();
}

/** @brief Measure the same coherent state used for spatial or spin-only VMC. */
template<typename TenElemT, typename QNT>
inline void RunPointGroupMeasure_(const heisenberg_params::PhysicalParams &physics,
                                 const qlpeps::MCMeasurementParams &measurement_params,
                                 const qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                 MPI_Comm comm) {
  using Sample = qlpeps::PointGroupProjectedSampleOBC<TenElemT, QNT>;
  using Updater = qlpeps::MCUpdateSquareNNPointGroupExchangeOBC;
  using Model = qlpeps::PointGroupSquareXXZModelOBC;
  using Measurer = qlpeps::MCPEPSMeasurer<TenElemT, QNT, Updater, Model,
                                        qlpeps::BMPSContractor, Sample>;
  const auto &projection = measurement_params.mc_params.point_group_projection;
  RequireProjectedConfiguration_(measurement_params.mc_params.initial_config, projection, comm);
  LogPointGroupProjection_(projection, comm);
  const Model model = MakePointGroupModel_(physics);
  Measurer measurer(sitps, measurement_params, comm, model, Updater{});
  measurer.Execute();
}

// ---------------- Measurement dispatcher (by model) ----------------
template<typename TenElemT, typename QNT>
inline void RunMeasureByModelOBC_(const heisenberg_params::PhysicalParams &phys,
                                  const heisenberg_params::AxisUpdateParams &axis_update,
                                  const heisenberg_params::MonteCarloNumericalParams &mc_params,
                                  const qlpeps::MCMeasurementParams &measurement_params,
                                  qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                                  MPI_Comm comm) {
  const std::string model_type = phys.ModelType.empty() ? "SquareHeisenberg" : phys.ModelType;
  const double j2 = phys.J2;
  // Same updater choice as VMC (VisitOBCUpdater).
  const auto run = [&](const auto &solver) {
    VisitOBCUpdater<QNT>(axis_update, mc_params, [&](auto updater) {
      heisenberg_vmcpeps::detail::ExecuteMeasure_<TenElemT, QNT, decltype(updater)>(
          sitps, measurement_params, comm, solver, std::move(updater));
    });
  };

  if (model_type == "SquareHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
      run(Model{});
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC;
      Model solver(j2);
      run(solver);
    }
    return;
  }

  if (model_type == "SquareXY") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::SquareSpinOneHalfXXZModelOBC; // XY: jz=0, jxy=1
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*pinning=*/0.0);
      run(solver);
    } else {
      using Model = qlpeps::SquareSpinOneHalfJ1J2XXZModelOBC; // XY with J2
      Model solver(/*jz=*/0.0, /*jxy=*/1.0, /*jz2=*/0.0, /*jxy2=*/j2, /*pinning=*/0.0);
      run(solver);
    }
    return;
  }

  if (model_type == "TriangleHeisenberg") {
    if (std::abs(j2) < 1e-15) {
      using Model = qlpeps::TriangularSpinOneHalfHeisenbergModelOBC;
      run(Model{});
    } else {
      using Model = qlpeps::TriangularSpinOneHalfJ1J2HeisenbergModelOBC;
      Model solver(j2);
      run(solver);
    }
    return;
  }

  // Fallback
  {
    using Model = qlpeps::SquareSpinOneHalfXXZModelOBC;
    run(Model{});
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
  if (params.point_group_projection.group != "None") {
    return heisenberg_vmcpeps::detail::RunPointGroupVmc_<TenElemT, QNT>(params, sitps, comm, rank);
  }
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

/**
 * @brief Run the measurement of the physics model; OBC sampling uses the updater VisitOBCUpdater
 *        selects from @p axis_update and @p mc_params (PBC ignores them: an enabled axis update
 *        is rejected for PBC when the parameters are parsed).
 */
template<typename TenElemT, typename QNT>
inline void RunMeasureByModel(const heisenberg_params::PhysicalParams &phys,
                              const heisenberg_params::AxisUpdateParams &axis_update,
                              const heisenberg_params::MonteCarloNumericalParams &mc_params,
                              const qlpeps::MCMeasurementParams &measurement_params,
                              qlpeps::SplitIndexTPS<TenElemT, QNT> &sitps,
                              MPI_Comm comm) {
  if (heisenberg_params::HasProjection(measurement_params.mc_params.point_group_projection)) {
    return heisenberg_vmcpeps::detail::RunPointGroupMeasure_<TenElemT, QNT>(
        phys, measurement_params, sitps, comm);
  }
  const bool is_pbc = (phys.BoundaryCondition == qlpeps::BoundaryCondition::Periodic);
  if (is_pbc) {
    heisenberg_vmcpeps::detail::RunMeasureByModelPBC_<TenElemT, QNT>(phys, measurement_params, sitps, comm);
  } else {
    heisenberg_vmcpeps::detail::RunMeasureByModelOBC_<TenElemT, QNT>(
        phys, axis_update, mc_params, measurement_params, sitps, comm);
  }
}

#endif // HEISENBERGVMCPEPS_MODEL_UPDATER_FACTORY_H
