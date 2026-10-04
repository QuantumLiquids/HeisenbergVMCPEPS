// SPDX-License-Identifier: LGPL-3.0-only
#ifndef HEISENBERGVMCPEPS_PROJECTION_PARAMS_H
#define HEISENBERGVMCPEPS_PROJECTION_PARAMS_H

#include "common_params.h"
#include "qlpeps/vmc_basic/point_group_projection.h"

namespace heisenberg_params {

/**
 * @brief Parse and validate coherent spatial projection and independent spin inversion.
 *
 * C4/D4 require a complete Lx=Ly square OBC lattice. Spin inversion additionally
 * requires a fixed Sz=0 sector. Group None retains the existing spin-only VMC
 * restrictions. C4 sectors 1 and 3 require a complex scalar build (RealCode=OFF).
 * D4 E is the central projector onto the two-dimensional isotypic subspace.
 */
inline qlpeps::PointGroupProjectionParams ReadProjectionParams(
    const char *algorithm_file, const PhysicalParams &physics,
    const MonteCarloNumericalParams &mc) {
  const auto values = qlpeps::config::ReadCaseParamsFile(algorithm_file);
  qlpeps::PointGroupProjectionParams params;
  params.group = qlpeps::config::ReadString(values, "PointGroup", "None");
  params.irrep = qlpeps::config::ReadString(
      values, "PointGroupIrrep", params.group == "C4" ? "0" : "A1");
  const double parity = qlpeps::config::ReadDouble(values, "SpinInversionParity", 0.0);
  if (parity != 0.0 && parity != 1.0 && parity != -1.0) {
    throw std::invalid_argument("SpinInversionParity must be 0, +1, or -1.");
  }
  params.spin_inversion_parity = static_cast<int>(parity);
  const bool c4 = params.group == "C4";
  const bool d4 = params.group == "D4";
  qlpeps::ValidatePointGroupProjectionParams(params);
#ifndef USE_COMPLEX
  if (c4 && (params.irrep == "1" || params.irrep == "3")) {
    throw std::invalid_argument("C4 sectors 1 and 3 require a complex build: configure with -DRealCode=OFF.");
  }
#endif
  if (!c4 && !d4 && params.spin_inversion_parity == 0) return params;
  if (physics.BoundaryCondition != qlpeps::BoundaryCondition::Open ||
      (physics.ModelType != "SquareHeisenberg" && physics.ModelType != "SquareXY") ||
      physics.RemoveCorner || physics.Lx < 2 || physics.Ly < 2) {
    throw std::invalid_argument("Projection requires a full square OBC lattice, Lx,Ly >= 2, and SquareHeisenberg or SquareXY.");
  }
  if ((c4 || d4) && !mc.MCRestrictU1) {
    throw std::invalid_argument("PointGroup projection currently requires MCRestrictU1=true: its exchange updater preserves spin counts.");
  }
  if ((c4 || d4) && physics.Lx != physics.Ly) {
    throw std::invalid_argument("C4/D4 point-group projection requires Lx == Ly.");
  }
  if (params.spin_inversion_parity != 0 &&
      (!mc.MCRestrictU1 || (physics.Lx * physics.Ly) % 2 != 0)) {
    throw std::invalid_argument("SpinInversionParity requires even Lx*Ly and MCRestrictU1=true (Sz=0).");
  }
  if (params.group == "None" && physics.J2 != 0.0) {
    throw std::invalid_argument("SpinInversionParity without PointGroup requires J2=0 for compatibility with spin-only VMC.");
  }
  return params;
}

/// Whether either coherent projector is enabled.
inline bool HasProjection(const qlpeps::PointGroupProjectionParams &params) {
  return params.group != "None" || params.spin_inversion_parity != 0;
}

/// Spatial projection requires the dedicated coherent-sum updater.
inline void ValidatePointGroupUpdater(const qlpeps::PointGroupProjectionParams &projection,
                                     const AxisUpdateParams &axis) {
  if (projection.group != "None" && axis.MCAxisUpdate) {
    throw std::invalid_argument("MCAxisUpdate=true is not supported with PointGroup projection.");
  }
}

}  // namespace heisenberg_params
#endif  // HEISENBERGVMCPEPS_PROJECTION_PARAMS_H
