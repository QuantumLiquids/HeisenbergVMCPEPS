// SPDX-License-Identifier: LGPL-3.0-only

/**
 * @file txt_to_tps_sitps.cpp
 * @brief Convert plain-text tensor blocks to TPS, then dump SplitIndexTPS.
 *
 * Input assumptions:
 * - One real number per non-empty line.
 * - Tensor blocks are separated by one or more blank lines.
 * - Each block is flattened from a source tensor index convention controlled by
 *   --src-index-order and --src-fast-order.
 *
 * Current default conversion (Yubin 2x2 unit cell):
 * - 4 tensor blocks.
 * - Source index labels: [phy, r, l, u, d].
 * - Source flatten order (fast -> slow): d, u, l, r, phy.
 * - Target TPS index order is fixed by qlpeps:
 *   (West, South, East, North, Physical) = (l, d, r, u, phy).
 */

#include "./qldouble.h"
#include "qlpeps/api/conversions.h"
#include "qlpeps/two_dim_tn/common/boundary_condition.h"
#include "qlpeps/two_dim_tn/tps/split_index_tps.h"
#include "qlpeps/two_dim_tn/tps/tps.h"

#include <array>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using BoundaryCondition = qlpeps::BoundaryCondition;
using TPS = qlpeps::TPS<TenElemT, QNT>;
using SITPS = qlpeps::SplitIndexTPS<TenElemT, QNT>;

enum Axis : size_t {
  kPhy = 0,
  kR = 1,
  kL = 2,
  kU = 3,
  kD = 4,
  kAxisCount = 5
};

struct CliOptions {
  std::string input_txt;
  std::string output_sitps_dir;
  size_t ly = 2;
  size_t lx = 2;
  size_t phy_dim = 2;
  size_t dim_r = 5;
  size_t dim_l = 5;
  size_t dim_u = 5;
  size_t dim_d = 5;
  std::string block_site_map = "0,0;0,1;1,0;1,1";
  std::string src_index_order = "phy,r,l,u,d";
  std::string src_fast_order = "d,u,l,r,phy";
  std::string boundary = "Periodic";
  std::optional<std::string> dump_tps_dir;
};

struct ParsedText {
  std::vector<std::vector<double>> blocks;
  size_t nonempty_lines = 0;
};

std::string Trim(const std::string &s) {
  size_t begin = 0;
  while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin]))) {
    ++begin;
  }
  size_t end = s.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(begin, end - begin);
}

std::string ToLower(std::string s) {
  for (char &ch : s) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }
  return s;
}

std::vector<std::string> Split(const std::string &s, const char delimiter) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delimiter)) {
    out.push_back(item);
  }
  return out;
}

bool ParsePositiveSizeT(const std::string &text, size_t *value) {
  if (value == nullptr) {
    return false;
  }
  if (text.empty()) {
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
    if (pos != text.size() || parsed == 0 || parsed > std::numeric_limits<size_t>::max()) {
      return false;
    }
    *value = static_cast<size_t>(parsed);
  } catch (...) {
    return false;
  }
  return true;
}

bool ParseNonNegativeSizeT(const std::string &text, size_t *value) {
  if (value == nullptr) {
    return false;
  }
  if (text.empty()) {
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
    if (pos != text.size() || parsed > std::numeric_limits<size_t>::max()) {
      return false;
    }
    *value = static_cast<size_t>(parsed);
  } catch (...) {
    return false;
  }
  return true;
}

bool ParseDoubleStrict(const std::string &text, double *value) {
  if (value == nullptr) {
    return false;
  }
  if (text.empty()) {
    return false;
  }
  errno = 0;
  char *end = nullptr;
  const char *begin = text.c_str();
  const double parsed = std::strtod(begin, &end);
  if (begin == end || errno == ERANGE) {
    return false;
  }
  while (*end != '\0') {
    if (!std::isspace(static_cast<unsigned char>(*end))) {
      return false;
    }
    ++end;
  }
  *value = parsed;
  return true;
}

bool CheckedMul(const size_t a, const size_t b, size_t *out) {
  if (out == nullptr) {
    return false;
  }
  if (a == 0 || b == 0) {
    *out = 0;
    return true;
  }
  if (a > std::numeric_limits<size_t>::max() / b) {
    return false;
  }
  *out = a * b;
  return true;
}

bool ParseAxisName(const std::string &name, Axis *axis) {
  if (axis == nullptr) {
    return false;
  }
  const std::string lowered = ToLower(Trim(name));
  if (lowered == "phy") {
    *axis = kPhy;
    return true;
  }
  if (lowered == "r") {
    *axis = kR;
    return true;
  }
  if (lowered == "l") {
    *axis = kL;
    return true;
  }
  if (lowered == "u") {
    *axis = kU;
    return true;
  }
  if (lowered == "d") {
    *axis = kD;
    return true;
  }
  return false;
}

bool ParseAxisOrder(
    const std::string &raw,
    const std::string &flag_name,
    std::vector<Axis> *order,
    std::string *err_msg) {
  if (order == nullptr || err_msg == nullptr) {
    return false;
  }
  const auto tokens = Split(raw, ',');
  if (tokens.size() != kAxisCount) {
    *err_msg = flag_name + " must have exactly 5 labels.";
    return false;
  }
  std::array<bool, kAxisCount> used{};
  order->clear();
  order->reserve(kAxisCount);
  for (const auto &token_raw : tokens) {
    Axis axis = kPhy;
    if (!ParseAxisName(token_raw, &axis)) {
      *err_msg = flag_name + " has unknown axis label: '" + Trim(token_raw) + "'.";
      return false;
    }
    if (used[axis]) {
      *err_msg = flag_name + " has duplicate axis label: '" + ToLower(Trim(token_raw)) + "'.";
      return false;
    }
    used[axis] = true;
    order->push_back(axis);
  }
  for (bool hit : used) {
    if (!hit) {
      *err_msg = flag_name + " is missing one or more axis labels.";
      return false;
    }
  }
  return true;
}

bool ParseBoundary(const std::string &raw, BoundaryCondition *bc, std::string *err_msg) {
  if (bc == nullptr || err_msg == nullptr) {
    return false;
  }
  const std::string lowered = ToLower(Trim(raw));
  if (lowered == "open" || lowered == "obc") {
    *bc = BoundaryCondition::Open;
    return true;
  }
  if (lowered == "periodic" || lowered == "pbc") {
    *bc = BoundaryCondition::Periodic;
    return true;
  }
  *err_msg = "--boundary must be Open/OBC or Periodic/PBC.";
  return false;
}

bool ParseBlockSiteMap(
    const std::string &raw,
    const size_t expected_count,
    const size_t ly,
    const size_t lx,
    std::vector<std::pair<size_t, size_t>> *site_map,
    std::string *err_msg) {
  if (site_map == nullptr || err_msg == nullptr) {
    return false;
  }
  site_map->clear();
  std::set<std::pair<size_t, size_t>> seen;
  const auto entries = Split(raw, ';');
  for (const auto &entry_raw : entries) {
    const std::string entry = Trim(entry_raw);
    if (entry.empty()) {
      continue;
    }
    const auto rc = Split(entry, ',');
    if (rc.size() != 2) {
      *err_msg = "--block-site-map entry must be 'row,col', got '" + entry + "'.";
      return false;
    }
    size_t row = 0;
    size_t col = 0;
    if (!ParseNonNegativeSizeT(Trim(rc[0]), &row) || !ParseNonNegativeSizeT(Trim(rc[1]), &col)) {
      *err_msg = "--block-site-map has invalid row/col in '" + entry + "'.";
      return false;
    }
    if (row >= ly || col >= lx) {
      *err_msg = "--block-site-map site out of range for lattice " + std::to_string(ly) + "x" +
                 std::to_string(lx) + ": '" + entry + "'.";
      return false;
    }
    const auto site = std::make_pair(row, col);
    if (!seen.insert(site).second) {
      *err_msg = "--block-site-map has duplicate site '" + entry + "'.";
      return false;
    }
    site_map->push_back(site);
  }
  if (site_map->size() != expected_count) {
    *err_msg = "--block-site-map count mismatch: got " + std::to_string(site_map->size()) +
               ", expected " + std::to_string(expected_count) + ".";
    return false;
  }
  return true;
}

void PrintUsage(const char *program) {
  std::cout
      << "Usage: " << program << " --input-txt <file> --output-sitps-dir <dir> [options]\n"
      << "Options:\n"
      << "  --ly <n>                  Lattice rows (default: 2)\n"
      << "  --lx <n>                  Lattice cols (default: 2)\n"
      << "  --phy-dim <n>             Physical dimension (default: 2)\n"
      << "  --dim-r <n>               Source r dimension (default: 5)\n"
      << "  --dim-l <n>               Source l dimension (default: 5)\n"
      << "  --dim-u <n>               Source u dimension (default: 5)\n"
      << "  --dim-d <n>               Source d dimension (default: 5)\n"
      << "  --block-site-map <spec>   Block->site map, 0-based row,col list\n"
      << "                            (default: \"0,0;0,1;1,0;1,1\")\n"
      << "  --src-index-order <spec>  Source axis labels, comma-separated\n"
      << "                            (default: \"phy,r,l,u,d\")\n"
      << "  --src-fast-order <spec>   Fast->slow flatten order\n"
      << "                            (default: \"d,u,l,r,phy\")\n"
      << "  --boundary <Open|Periodic> Boundary condition (default: Periodic)\n"
      << "  --dump-tps-dir <dir>      Optional intermediate TPS dump dir\n"
      << "  --help                    Show this message\n";
}

bool ParseCli(int argc, char **argv, CliOptions *opts, std::string *err_msg) {
  if (opts == nullptr || err_msg == nullptr) {
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
    const std::string value = argv[++i];
    if (arg == "--input-txt") {
      opts->input_txt = value;
    } else if (arg == "--output-sitps-dir") {
      opts->output_sitps_dir = value;
    } else if (arg == "--ly") {
      if (!ParsePositiveSizeT(value, &opts->ly)) {
        *err_msg = "Invalid --ly: " + value;
        return false;
      }
    } else if (arg == "--lx") {
      if (!ParsePositiveSizeT(value, &opts->lx)) {
        *err_msg = "Invalid --lx: " + value;
        return false;
      }
    } else if (arg == "--phy-dim") {
      if (!ParsePositiveSizeT(value, &opts->phy_dim)) {
        *err_msg = "Invalid --phy-dim: " + value;
        return false;
      }
    } else if (arg == "--dim-r") {
      if (!ParsePositiveSizeT(value, &opts->dim_r)) {
        *err_msg = "Invalid --dim-r: " + value;
        return false;
      }
    } else if (arg == "--dim-l") {
      if (!ParsePositiveSizeT(value, &opts->dim_l)) {
        *err_msg = "Invalid --dim-l: " + value;
        return false;
      }
    } else if (arg == "--dim-u") {
      if (!ParsePositiveSizeT(value, &opts->dim_u)) {
        *err_msg = "Invalid --dim-u: " + value;
        return false;
      }
    } else if (arg == "--dim-d") {
      if (!ParsePositiveSizeT(value, &opts->dim_d)) {
        *err_msg = "Invalid --dim-d: " + value;
        return false;
      }
    } else if (arg == "--block-site-map") {
      opts->block_site_map = value;
    } else if (arg == "--src-index-order") {
      opts->src_index_order = value;
    } else if (arg == "--src-fast-order") {
      opts->src_fast_order = value;
    } else if (arg == "--boundary") {
      opts->boundary = value;
    } else if (arg == "--dump-tps-dir") {
      opts->dump_tps_dir = value;
    } else {
      *err_msg = "Unknown argument: " + arg;
      return false;
    }
  }
  if (opts->input_txt.empty()) {
    *err_msg = "--input-txt is required.";
    return false;
  }
  if (opts->output_sitps_dir.empty()) {
    *err_msg = "--output-sitps-dir is required.";
    return false;
  }
  return true;
}

bool ParseTxtBlocks(const std::string &path, ParsedText *parsed, std::string *err_msg) {
  if (parsed == nullptr || err_msg == nullptr) {
    return false;
  }
  parsed->blocks.clear();
  parsed->nonempty_lines = 0;

  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    *err_msg = "Failed to open input TXT: " + path;
    return false;
  }

  std::vector<double> current_block;
  std::string line;
  size_t line_number = 0;
  while (std::getline(ifs, line)) {
    ++line_number;
    const std::string trimmed = Trim(line);
    if (trimmed.empty()) {
      if (!current_block.empty()) {
        parsed->blocks.push_back(current_block);
        current_block.clear();
      }
      continue;
    }
    double value = 0.0;
    if (!ParseDoubleStrict(trimmed, &value)) {
      *err_msg = "Malformed real number at line " + std::to_string(line_number) + ": '" + trimmed + "'.";
      return false;
    }
    ++parsed->nonempty_lines;
    current_block.push_back(value);
  }
  if (!current_block.empty()) {
    parsed->blocks.push_back(current_block);
  }
  if (parsed->blocks.empty()) {
    *err_msg = "Input TXT contains no tensor data blocks.";
    return false;
  }
  return true;
}

std::array<size_t, kAxisCount> DecodeFlatIndex(
    const size_t flat,
    const std::vector<Axis> &fast_order,
    const std::array<size_t, kAxisCount> &dims_by_axis) {
  std::array<size_t, kAxisCount> coord{};
  size_t remain = flat;
  for (Axis axis : fast_order) {
    const size_t dim = dims_by_axis[axis];
    coord[axis] = remain % dim;
    remain /= dim;
  }
  return coord;
}

}  // namespace

int main(int argc, char **argv) {
  CliOptions opts;
  std::string err_msg;
  if (!ParseCli(argc, argv, &opts, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    PrintUsage(argv[0]);
    return -1;
  }

  BoundaryCondition bc = BoundaryCondition::Periodic;
  if (!ParseBoundary(opts.boundary, &bc, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    return -2;
  }

  std::vector<Axis> src_index_order;
  if (!ParseAxisOrder(opts.src_index_order, "--src-index-order", &src_index_order, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    return -3;
  }
  std::vector<Axis> src_fast_order;
  if (!ParseAxisOrder(opts.src_fast_order, "--src-fast-order", &src_fast_order, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    return -4;
  }

  size_t expected_blocks = 0;
  if (!CheckedMul(opts.ly, opts.lx, &expected_blocks)) {
    std::cerr << "ERROR: ly * lx overflow." << std::endl;
    return -5;
  }

  std::vector<std::pair<size_t, size_t>> block_site_map;
  if (!ParseBlockSiteMap(
          opts.block_site_map, expected_blocks, opts.ly, opts.lx, &block_site_map, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    return -6;
  }

  ParsedText parsed;
  if (!ParseTxtBlocks(opts.input_txt, &parsed, &err_msg)) {
    std::cerr << "ERROR: " << err_msg << std::endl;
    return -7;
  }
  if (parsed.blocks.size() != expected_blocks) {
    std::cerr << "ERROR: Block count mismatch. Parsed " << parsed.blocks.size()
              << " blocks, expected " << expected_blocks << " (ly*lx)." << std::endl;
    return -8;
  }

  const std::array<size_t, kAxisCount> dims_by_axis = {
      opts.phy_dim, opts.dim_r, opts.dim_l, opts.dim_u, opts.dim_d};
  size_t expected_values_per_block = 1;
  for (Axis axis : src_index_order) {
    if (!CheckedMul(expected_values_per_block, dims_by_axis[axis], &expected_values_per_block)) {
      std::cerr << "ERROR: Overflow while computing expected values per block." << std::endl;
      return -9;
    }
  }

  for (size_t i = 0; i < parsed.blocks.size(); ++i) {
    if (parsed.blocks[i].size() != expected_values_per_block) {
      std::cerr << "ERROR: Block " << (i + 1) << " value count mismatch. Got "
                << parsed.blocks[i].size() << ", expected " << expected_values_per_block << "."
                << std::endl;
      return -10;
    }
  }

  std::cout << "Parsed TXT blocks: " << parsed.blocks.size()
            << ", values per block: " << expected_values_per_block << std::endl;

  qlten::hp_numeric::SetTensorManipulationThreads(1);

  const IndexT west_index({QNSctT(qn0, opts.dim_l)}, qlten::TenIndexDirType::IN);
  const IndexT south_index({QNSctT(qn0, opts.dim_d)}, qlten::TenIndexDirType::OUT);
  const IndexT east_index({QNSctT(qn0, opts.dim_r)}, qlten::TenIndexDirType::OUT);
  const IndexT north_index({QNSctT(qn0, opts.dim_u)}, qlten::TenIndexDirType::IN);

#ifdef U1SYM
  if (opts.phy_dim != pb_out.dim()) {
    std::cerr << "ERROR: U1 build currently supports only --phy-dim=" << pb_out.dim() << "." << std::endl;
    return -11;
  }
  const IndexT physical_index = pb_out;
#else
  const IndexT physical_index({QNSctT(qn0, opts.phy_dim)}, qlten::TenIndexDirType::OUT);
#endif

  TPS tps(opts.ly, opts.lx, bc);
  for (size_t block = 0; block < parsed.blocks.size(); ++block) {
    const auto [row, col] = block_site_map[block];
    Tensor local_ten({west_index, south_index, east_index, north_index, physical_index});
    const auto &values = parsed.blocks[block];
    for (size_t flat = 0; flat < values.size(); ++flat) {
      const auto coord = DecodeFlatIndex(flat, src_fast_order, dims_by_axis);
      local_ten({coord[kL], coord[kD], coord[kR], coord[kU], coord[kPhy]}) = values[flat];
    }
    tps({row, col}) = local_ten;
  }

  if (opts.dump_tps_dir.has_value()) {
    tps.Dump(*opts.dump_tps_dir);
    std::cout << "Intermediate TPS dumped to: " << *opts.dump_tps_dir << std::endl;
  }

  SITPS sitps = qlpeps::ToSplitIndexTPS<TenElemT, QNT>(tps);
  sitps.Dump(opts.output_sitps_dir);
  std::cout << "SplitIndexTPS dumped to: " << opts.output_sitps_dir << std::endl;

  SITPS loadback;
  if (!loadback.Load(opts.output_sitps_dir)) {
    std::cerr << "ERROR: Loadback check failed for dumped SITPS: " << opts.output_sitps_dir << std::endl;
    return -12;
  }
  if (loadback.rows() != opts.ly || loadback.cols() != opts.lx) {
    std::cerr << "ERROR: Loadback lattice shape mismatch. Got " << loadback.rows() << "x" << loadback.cols()
              << ", expected " << opts.ly << "x" << opts.lx << "." << std::endl;
    return -13;
  }
  if (loadback.PhysicalDim() != opts.phy_dim) {
    std::cerr << "ERROR: Loadback physical dimension mismatch. Got " << loadback.PhysicalDim()
              << ", expected " << opts.phy_dim << "." << std::endl;
    return -14;
  }
  if (loadback.GetBoundaryCondition() != bc) {
    std::cerr << "ERROR: Loadback boundary mismatch." << std::endl;
    return -15;
  }
  std::cout << "Loadback check passed." << std::endl;
  return 0;
}
