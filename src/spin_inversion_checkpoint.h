// SPDX-License-Identifier: LGPL-3.0-only
#ifndef SPIN_INVERSION_CHECKPOINT_H
#define SPIN_INVERSION_CHECKPOINT_H

#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>
#include "qlpeps/algorithm/vmc_update/vmc_peps_optimizer_params.h"

namespace spin_inversion_io {

/// Read projected-state metadata; unmarked legacy PEPS directories return zero.
inline int ReadParity(const std::filesystem::path &directory) {
  auto normalized = directory.lexically_normal();
  if (normalized.filename().empty()) normalized = normalized.parent_path();
  auto marker = normalized / "spin_inversion_parity.txt";
  // Periodic step directories inherit the marker written before optimization.
  if (!std::filesystem::exists(marker) &&
      normalized.filename().string().starts_with("step_")) {
    marker = normalized.parent_path() / "spin_inversion_parity.txt";
  }
  if (!std::filesystem::exists(marker)) return 0;
  std::ifstream input(marker);
  int parity = 0;
  std::string extra;
  if (!(input >> parity) || (parity != 1 && parity != -1) || (input >> extra)) {
    throw std::runtime_error("Invalid spin-inversion metadata: " + marker.string());
  }
  return parity;
}

/// Refuse to reinterpret marked base tensors as a different physical ansatz.
inline void RequireParity(const std::string &directory, int requested_parity) {
  const int saved_parity = ReadParity(directory);
  if (saved_parity != 0 && saved_parity != requested_parity) {
    throw std::runtime_error(
        "SpinInversionParity mismatch for '" + directory + "': saved " +
        std::to_string(saved_parity) + ", requested " +
        std::to_string(requested_parity) +
        ". These are base PEPS tensors of a projected wavefunction; plain "
        "mc_measure does not support this state. Preserve its metadata.");
  }
}

/// Run filesystem operations on rank zero and propagate any error to all ranks.
inline void Collective(MPI_Comm comm, int rank, const std::function<void()> &operation) {
  std::string error;
  if (rank == 0) {
    try {
      operation();
    } catch (const std::exception &exception) {
      error = exception.what();
    }
  }
  int length = static_cast<int>(error.size());
  ::MPI_Bcast(&length, 1, MPI_INT, 0, comm);
  error.resize(static_cast<size_t>(length));
  if (length != 0) {
    ::MPI_Bcast(error.data(), length, MPI_CHAR, 0, comm);
    throw std::runtime_error(error);
  }
}

/// Mark final/lowest and periodic output roots before any projected state is dumped.
inline void PrepareOutputs(const qlpeps::VMCPEPSOptimizerParams &params, int parity,
                           MPI_Comm comm, int rank) {
  Collective(comm, rank, [&] {
    std::vector<std::string> directories;
    if (!params.tps_dump_base_name.empty()) {
      directories.push_back(params.tps_dump_base_name + "final");
      directories.push_back(params.tps_dump_base_name + "lowest");
    }
    const auto &checkpoint = params.optimizer_params.checkpoint_params;
    if (checkpoint.IsEnabled()) directories.push_back(checkpoint.base_path);

    // Validate every destination before creating or marking any output directory.
    for (const auto &directory : directories) RequireParity(directory, parity);
    if (parity != 0 && checkpoint.IsEnabled() &&
        ReadParity(checkpoint.base_path) == 0 &&
        std::filesystem::exists(checkpoint.base_path)) {
      for (const auto &entry : std::filesystem::directory_iterator(checkpoint.base_path)) {
        if (entry.is_directory() && entry.path().filename().string().starts_with("step_")) {
          throw std::runtime_error(
              "Cannot mark checkpoint root '" + checkpoint.base_path +
              "' as projected: it contains existing unmarked step_* directories. "
              "Choose a new CheckpointBasePath to preserve the old checkpoints' meaning.");
        }
      }
    }
    if (parity == 0) return;
    for (const auto &directory : directories) {
      std::filesystem::create_directories(directory);
      std::ofstream output(std::filesystem::path(directory) / "spin_inversion_parity.txt");
      output << parity << '\n';
      output.close();
      if (!output) throw std::runtime_error("Cannot write spin-inversion metadata in " + directory);
    }
  });
}

}  // namespace spin_inversion_io
#endif  // SPIN_INVERSION_CHECKPOINT_H
