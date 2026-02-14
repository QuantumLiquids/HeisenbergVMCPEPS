// SPDX-License-Identifier: LGPL-3.0-only

/*
 * Author: Codex
 * Description: Tile/replicate SplitIndexTPS with strict leg-consistency checks.
 *
 * Usage:
 *   ./sitps_tile --input-dir <dir> --output-dir <dir>
 *                --target-ly <Ly> --target-lx <Lx>
 *                [--unit-ly <ly>] [--unit-lx <lx>]
 *                [--physics-json <file>]
 */

#include "./qldouble.h"
#include "common_params.h"
#include "qlpeps/two_dim_tn/tps/split_index_tps.h"

#include <array>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using SITPS = qlpeps::SplitIndexTPS<TenElemT, QNT>;
using BoundaryCondition = qlpeps::BoundaryCondition;

struct CliOptions {
  std::string input_dir;
  std::string output_dir;
  size_t target_ly = 0;
  size_t target_lx = 0;
  std::optional<size_t> unit_ly;
  std::optional<size_t> unit_lx;
  std::optional<std::string> physics_json;
};

class ValidationReporter {
 public:
  explicit ValidationReporter(size_t max_errors) : max_errors_(max_errors) {}

  void Add(const std::string &message) {
    if (errors_.size() < max_errors_) {
      errors_.push_back(message);
    } else {
      truncated_ = true;
    }
  }

  bool HasErrors() const { return !errors_.empty() || truncated_; }

  void Print(const std::string &title) const {
    std::cerr << "ERROR: " << title << std::endl;
    for (const auto &err : errors_) {
      std::cerr << "  - " << err << std::endl;
    }
    if (truncated_) {
      std::cerr << "  - ... more mismatches omitted ..." << std::endl;
    }
  }

 private:
  size_t max_errors_;
  bool truncated_ = false;
  std::vector<std::string> errors_;
};

bool ParseSizeT(const std::string &text, size_t *value) {
  if (value == nullptr || text.empty()) {
    return false;
  }
  for (char ch : text) {
    if (!std::isdigit(static_cast<unsigned char>(ch))) {
      return false;
    }
  }
  size_t pos = 0;
  try {
    const unsigned long long parsed = std::stoull(text, &pos, 10);
    if (pos != text.size()) {
      return false;
    }
    if (parsed > std::numeric_limits<size_t>::max()) {
      return false;
    }
    *value = static_cast<size_t>(parsed);
  } catch (...) {
    return false;
  }
  return true;
}

std::string BoundaryConditionToString(const BoundaryCondition bc) {
  return (bc == BoundaryCondition::Periodic) ? "Periodic" : "Open";
}

std::string IndexToString(const IndexT &idx) {
  std::ostringstream oss;
  oss << "{dim=" << idx.dim() << ", dir=";
  const auto dir = idx.GetDir();
  if (dir == qlten::TenIndexDirType::IN) {
    oss << "IN";
  } else if (dir == qlten::TenIndexDirType::OUT) {
    oss << "OUT";
  } else {
    oss << "NDIR";
  }
  oss << ", qn=[";
  for (size_t i = 0; i < idx.GetQNSctNum(); ++i) {
    if (i > 0) {
      oss << ", ";
    }
    const auto &sct = idx.GetQNSct(i);
    oss << sct.GetQn() << ":" << sct.dim();
  }
  oss << "]}";
  return oss.str();
}

const Tensor *GetRepresentativeTensor(const std::vector<Tensor> &site_comps) {
  for (const auto &ten : site_comps) {
    if (!ten.IsDefault() && ten.GetIndexes().size() >= 4) {
      return &ten;
    }
  }
  return nullptr;
}

void PrintUsage(const char *prog) {
  std::cout
      << "Usage: " << prog << " --input-dir <dir> --output-dir <dir>"
      << " --target-ly <Ly> --target-lx <Lx>"
      << " [--unit-ly <ly>] [--unit-lx <lx>] [--physics-json <file>]\n";
}

bool ParseCli(int argc, char **argv, CliOptions *opts, std::string *err_msg) {
  if (opts == nullptr || err_msg == nullptr) {
    return false;
  }
  if (argc == 1) {
    *err_msg = "No arguments provided.";
    return false;
  }
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      std::exit(0);
    }
    if (i + 1 >= argc) {
      *err_msg = "Missing value for argument: " + arg;
      return false;
    }
    const std::string val = argv[++i];
    if (arg == "--input-dir") {
      opts->input_dir = val;
    } else if (arg == "--output-dir") {
      opts->output_dir = val;
    } else if (arg == "--target-ly") {
      if (!ParseSizeT(val, &opts->target_ly)) {
        *err_msg = "Invalid --target-ly: " + val;
        return false;
      }
    } else if (arg == "--target-lx") {
      if (!ParseSizeT(val, &opts->target_lx)) {
        *err_msg = "Invalid --target-lx: " + val;
        return false;
      }
    } else if (arg == "--unit-ly") {
      size_t parsed = 0;
      if (!ParseSizeT(val, &parsed)) {
        *err_msg = "Invalid --unit-ly: " + val;
        return false;
      }
      opts->unit_ly = parsed;
    } else if (arg == "--unit-lx") {
      size_t parsed = 0;
      if (!ParseSizeT(val, &parsed)) {
        *err_msg = "Invalid --unit-lx: " + val;
        return false;
      }
      opts->unit_lx = parsed;
    } else if (arg == "--physics-json") {
      opts->physics_json = val;
    } else {
      *err_msg = "Unknown argument: " + arg;
      return false;
    }
  }
  if (opts->input_dir.empty() || opts->output_dir.empty()) {
    *err_msg = "Both --input-dir and --output-dir are required.";
    return false;
  }
  if (opts->target_ly == 0 || opts->target_lx == 0) {
    *err_msg = "--target-ly and --target-lx must be positive.";
    return false;
  }
  return true;
}

bool IsTrgSupportedLinearSize(size_t n) {
  if (n == 0) {
    return false;
  }
  while ((n % 2ULL) == 0ULL) {
    n /= 2ULL;
  }
  return (n == 1ULL) || (n == 3ULL);
}

bool ValidateSiteLevel(const SITPS &sitps, ValidationReporter *reporter) {
  if (reporter == nullptr) {
    return false;
  }
  size_t expected_phy_dim = 0;
  for (size_t row = 0; row < sitps.rows(); ++row) {
    for (size_t col = 0; col < sitps.cols(); ++col) {
      const auto &site = sitps({row, col});
      if (site.empty()) {
        reporter->Add("site(" + std::to_string(row) + "," + std::to_string(col) + ") has empty component list.");
        continue;
      }
      if (expected_phy_dim == 0) {
        expected_phy_dim = site.size();
      } else if (site.size() != expected_phy_dim) {
        reporter->Add(
            "site(" + std::to_string(row) + "," + std::to_string(col) + ") physical dimension mismatch: got " +
            std::to_string(site.size()) + ", expected " + std::to_string(expected_phy_dim) + ".");
      }

      bool has_reference = false;
      size_t ref_rank = 0;
      std::array<IndexT, 4> ref_virtual{};
      for (size_t comp = 0; comp < site.size(); ++comp) {
        const Tensor &ten = site[comp];
        if (ten.IsDefault()) {
          reporter->Add("site(" + std::to_string(row) + "," + std::to_string(col) + "), comp " +
                        std::to_string(comp) + " is default/uninitialized.");
          continue;
        }
        const auto &idxes = ten.GetIndexes();
        if (idxes.size() < 4) {
          reporter->Add("site(" + std::to_string(row) + "," + std::to_string(col) + "), comp " +
                        std::to_string(comp) + " rank < 4.");
          continue;
        }
        if (!has_reference) {
          ref_rank = idxes.size();
          for (size_t leg = 0; leg < 4; ++leg) {
            ref_virtual[leg] = idxes[leg];
          }
          has_reference = true;
        } else {
          if (idxes.size() != ref_rank) {
            reporter->Add("site(" + std::to_string(row) + "," + std::to_string(col) + "), comp " +
                          std::to_string(comp) + " rank mismatch: got " + std::to_string(idxes.size()) +
                          ", expected " + std::to_string(ref_rank) + ".");
          }
          for (size_t leg = 0; leg < 4; ++leg) {
            if (idxes[leg] != ref_virtual[leg]) {
              reporter->Add("site(" + std::to_string(row) + "," + std::to_string(col) + "), comp " +
                            std::to_string(comp) + ", leg " + std::to_string(leg) +
                            " differs from site reference.");
            }
          }
        }
      }
    }
  }
  return !reporter->HasErrors();
}

bool ValidateBondLevel(const SITPS &sitps, ValidationReporter *reporter) {
  if (reporter == nullptr) {
    return false;
  }
  const bool is_pbc = (sitps.GetBoundaryCondition() == BoundaryCondition::Periodic);

  auto check_horizontal = [&](size_t row, size_t col_a, size_t col_b) {
    const auto &site_a = sitps({row, col_a});
    const auto &site_b = sitps({row, col_b});
    const Tensor *ten_a = GetRepresentativeTensor(site_a);
    const Tensor *ten_b = GetRepresentativeTensor(site_b);
    if (ten_a == nullptr || ten_b == nullptr) {
      reporter->Add("cannot verify horizontal bond at row=" + std::to_string(row) + ", cols=(" +
                    std::to_string(col_a) + "," + std::to_string(col_b) +
                    "): representative tensor missing.");
      return;
    }
    const IndexT expected = qlten::InverseIndex(ten_b->GetIndex(0));
    const IndexT actual = ten_a->GetIndex(2);
    if (actual != expected) {
      reporter->Add("horizontal mismatch: A(" + std::to_string(row) + "," + std::to_string(col_a) + ").E != "
                    "Inverse(B(" + std::to_string(row) + "," + std::to_string(col_b) + ").W); actual " +
                    IndexToString(actual) + ", expected " + IndexToString(expected));
    }
  };

  auto check_vertical = [&](size_t row_a, size_t row_b, size_t col) {
    const auto &site_a = sitps({row_a, col});
    const auto &site_b = sitps({row_b, col});
    const Tensor *ten_a = GetRepresentativeTensor(site_a);
    const Tensor *ten_b = GetRepresentativeTensor(site_b);
    if (ten_a == nullptr || ten_b == nullptr) {
      reporter->Add("cannot verify vertical bond at col=" + std::to_string(col) + ", rows=(" +
                    std::to_string(row_a) + "," + std::to_string(row_b) +
                    "): representative tensor missing.");
      return;
    }
    const IndexT expected = qlten::InverseIndex(ten_b->GetIndex(3));
    const IndexT actual = ten_a->GetIndex(1);
    if (actual != expected) {
      reporter->Add("vertical mismatch: A(" + std::to_string(row_a) + "," + std::to_string(col) + ").S != "
                    "Inverse(B(" + std::to_string(row_b) + "," + std::to_string(col) + ").N); actual " +
                    IndexToString(actual) + ", expected " + IndexToString(expected));
    }
  };

  if (is_pbc) {
    for (size_t row = 0; row < sitps.rows(); ++row) {
      for (size_t col = 0; col < sitps.cols(); ++col) {
        check_horizontal(row, col, (col + 1) % sitps.cols());
        check_vertical(row, (row + 1) % sitps.rows(), col);
      }
    }
  } else {
    for (size_t row = 0; row < sitps.rows(); ++row) {
      for (size_t col = 0; col + 1 < sitps.cols(); ++col) {
        check_horizontal(row, col, col + 1);
      }
    }
    for (size_t row = 0; row + 1 < sitps.rows(); ++row) {
      for (size_t col = 0; col < sitps.cols(); ++col) {
        check_vertical(row, row + 1, col);
      }
    }
  }
  return !reporter->HasErrors();
}

bool ValidateOBCBoundary(const SITPS &sitps, ValidationReporter *reporter) {
  if (reporter == nullptr) {
    return false;
  }
  if (sitps.GetBoundaryCondition() != BoundaryCondition::Open) {
    return true;
  }

  for (size_t row = 0; row < sitps.rows(); ++row) {
    for (size_t col = 0; col < sitps.cols(); ++col) {
      const Tensor *ten = GetRepresentativeTensor(sitps({row, col}));
      if (ten == nullptr) {
        reporter->Add("cannot verify OBC boundary at site(" + std::to_string(row) + "," + std::to_string(col) +
                      "): representative tensor missing.");
        continue;
      }
      if (col == 0 && ten->GetIndex(0).dim() != 1) {
        reporter->Add("OBC boundary mismatch at site(" + std::to_string(row) + "," + std::to_string(col) +
                      "): left leg dim != 1.");
      }
      if (row + 1 == sitps.rows() && ten->GetIndex(1).dim() != 1) {
        reporter->Add("OBC boundary mismatch at site(" + std::to_string(row) + "," + std::to_string(col) +
                      "): bottom leg dim != 1.");
      }
      if (col + 1 == sitps.cols() && ten->GetIndex(2).dim() != 1) {
        reporter->Add("OBC boundary mismatch at site(" + std::to_string(row) + "," + std::to_string(col) +
                      "): right leg dim != 1.");
      }
      if (row == 0 && ten->GetIndex(3).dim() != 1) {
        reporter->Add("OBC boundary mismatch at site(" + std::to_string(row) + "," + std::to_string(col) +
                      "): top leg dim != 1.");
      }
    }
  }
  return !reporter->HasErrors();
}

bool ValidateSITPS(const SITPS &sitps, const std::string &name) {
  ValidationReporter reporter(24);
  ValidateSiteLevel(sitps, &reporter);
  ValidateBondLevel(sitps, &reporter);
  ValidateOBCBoundary(sitps, &reporter);
  if (reporter.HasErrors()) {
    reporter.Print("SITPS validation failed for '" + name + "'");
    return false;
  }
  return true;
}

std::pair<size_t, size_t> MapPBCSite(size_t dst_row, size_t dst_col, size_t unit_ly, size_t unit_lx) {
  return {dst_row % unit_ly, dst_col % unit_lx};
}

std::pair<size_t, size_t> MapOBCSite(
    size_t dst_row, size_t dst_col, size_t target_ly, size_t target_lx, size_t src_ly, size_t src_lx) {
  size_t src_row = 0;
  if (dst_row == 0) {
    src_row = 0;
  } else if (dst_row + 1 == target_ly) {
    src_row = src_ly - 1;
  } else {
    src_row = 1 + ((dst_row - 1) % (src_ly - 2));
  }

  size_t src_col = 0;
  if (dst_col == 0) {
    src_col = 0;
  } else if (dst_col + 1 == target_lx) {
    src_col = src_lx - 1;
  } else {
    src_col = 1 + ((dst_col - 1) % (src_lx - 2));
  }
  return {src_row, src_col};
}

bool ValidatePhysicsCrossCheck(const CliOptions &opts, BoundaryCondition bc_from_sitps) {
  if (!opts.physics_json.has_value()) {
    return true;
  }
  try {
    const heisenberg_params::PhysicalParams phys(opts.physics_json->c_str());
    if (phys.Ly != opts.target_ly) {
      std::cerr << "ERROR: physics Ly mismatch: physics=" << phys.Ly
                << ", target_ly=" << opts.target_ly << std::endl;
      return false;
    }
    if (phys.Lx != opts.target_lx) {
      std::cerr << "ERROR: physics Lx mismatch: physics=" << phys.Lx
                << ", target_lx=" << opts.target_lx << std::endl;
      return false;
    }
    if (phys.BoundaryCondition != bc_from_sitps) {
      std::cerr << "ERROR: physics BoundaryCondition mismatch: physics="
                << BoundaryConditionToString(phys.BoundaryCondition) << ", sitps="
                << BoundaryConditionToString(bc_from_sitps) << std::endl;
      return false;
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: Failed to parse --physics-json '" << *opts.physics_json
              << "': " << e.what() << std::endl;
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  CliOptions opts;
  std::string parse_err;
  if (!ParseCli(argc, argv, &opts, &parse_err)) {
    std::cerr << "ERROR: " << parse_err << std::endl;
    PrintUsage(argv[0]);
    return 1;
  }

  SITPS source;
  try {
    if (!source.Load(opts.input_dir)) {
      std::cerr << "ERROR: Failed to load SplitIndexTPS from '" << opts.input_dir << "'." << std::endl;
      return 2;
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: Exception while loading SplitIndexTPS from '" << opts.input_dir
              << "': " << e.what() << std::endl;
    return 2;
  }

  const size_t src_ly = source.rows();
  const size_t src_lx = source.cols();
  const BoundaryCondition bc = source.GetBoundaryCondition();

  if (opts.target_ly < src_ly || opts.target_lx < src_lx) {
    std::cerr << "ERROR: target size must be >= source size. source=(" << src_ly << "," << src_lx
              << "), target=(" << opts.target_ly << "," << opts.target_lx << ")." << std::endl;
    return 3;
  }

  if (!ValidatePhysicsCrossCheck(opts, bc)) {
    return 4;
  }

  if (!ValidateSITPS(source, "source")) {
    return 5;
  }

  size_t unit_ly = src_ly;
  size_t unit_lx = src_lx;
  if (bc == BoundaryCondition::Periodic) {
    unit_ly = opts.unit_ly.value_or(src_ly);
    unit_lx = opts.unit_lx.value_or(src_lx);
    if (unit_ly == 0 || unit_lx == 0 || unit_ly > src_ly || unit_lx > src_lx) {
      std::cerr << "ERROR: PBC unit cell must satisfy 1 <= unit <= source size. "
                << "source=(" << src_ly << "," << src_lx << "), unit=("
                << unit_ly << "," << unit_lx << ")." << std::endl;
      return 6;
    }
    if ((opts.target_ly % unit_ly) != 0 || (opts.target_lx % unit_lx) != 0) {
      std::cerr << "ERROR: PBC requires target dimensions to be integer multiples of unit cell. "
                << "target=(" << opts.target_ly << "," << opts.target_lx << "), unit=("
                << unit_ly << "," << unit_lx << ")." << std::endl;
      return 7;
    }
  } else {
    if (opts.unit_ly.has_value() || opts.unit_lx.has_value()) {
      std::cerr << "ERROR: --unit-ly/--unit-lx are not allowed for OBC mode." << std::endl;
      return 8;
    }
    if (src_ly < 3 || src_lx < 3) {
      std::cerr << "ERROR: OBC boundary-preserving tiling requires source size >= 3x3, got "
                << src_ly << "x" << src_lx << "." << std::endl;
      return 9;
    }
  }

  SITPS tiled;
  try {
    tiled = SITPS(opts.target_ly, opts.target_lx, bc);
    for (size_t row = 0; row < opts.target_ly; ++row) {
      for (size_t col = 0; col < opts.target_lx; ++col) {
        std::pair<size_t, size_t> src_site;
        if (bc == BoundaryCondition::Periodic) {
          src_site = MapPBCSite(row, col, unit_ly, unit_lx);
        } else {
          src_site = MapOBCSite(row, col, opts.target_ly, opts.target_lx, src_ly, src_lx);
        }
        tiled({row, col}) = source({src_site.first, src_site.second});
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: Failed to build tiled SITPS for target=(" << opts.target_ly
              << "," << opts.target_lx << "): " << e.what() << std::endl;
    return 12;
  }

  if (!ValidateSITPS(tiled, "tiled")) {
    return 10;
  }

  if (bc == BoundaryCondition::Periodic) {
    if (opts.target_ly != opts.target_lx) {
      std::cerr << "WARNING: PBC target is not square (" << opts.target_ly << "x" << opts.target_lx
                << "); TRG backend may reject it." << std::endl;
    }
    if (!IsTrgSupportedLinearSize(opts.target_ly) || !IsTrgSupportedLinearSize(opts.target_lx)) {
      std::cerr << "WARNING: PBC target size is not in L=2^k or L=3*2^k family; TRG backend may reject it."
                << std::endl;
    }
  }

  try {
    tiled.Dump(opts.output_dir);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: Failed to dump tiled SITPS to '" << opts.output_dir
              << "': " << e.what() << std::endl;
    return 11;
  }

  std::cout << "SITPS tiled successfully.\n"
            << "  Input:  " << opts.input_dir << " (" << src_ly << "x" << src_lx << ", "
            << BoundaryConditionToString(bc) << ")\n"
            << "  Output: " << opts.output_dir << " (" << opts.target_ly << "x" << opts.target_lx << ")\n";
  if (bc == BoundaryCondition::Periodic) {
    std::cout << "  Unit:   " << unit_ly << "x" << unit_lx << std::endl;
  }
  return 0;
}
