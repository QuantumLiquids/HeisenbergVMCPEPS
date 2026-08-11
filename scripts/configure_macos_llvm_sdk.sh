#!/usr/bin/env bash
set -euo pipefail

# Configure this project on macOS using Homebrew LLVM + OpenMP,
# while forcing the macOS SDK libc++ headers to avoid ABI mismatches
# (e.g. undefined symbol: std::__1::__hash_memory).

build_dir="${1:-build}"
build_type="${2:-Release}"

sdk="$(xcrun --sdk macosx --show-sdk-path)"

rm -rf "${build_dir}/CMakeCache.txt" "${build_dir}/CMakeFiles"

cmake -S . -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE="${build_type}" \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
  -DCMAKE_OSX_SYSROOT="${sdk}" \
  -DCMAKE_CXX_FLAGS="-nostdinc++ -isystem ${sdk}/usr/include/c++/v1"

echo "[ok] Configured: ${build_dir} (CMAKE_BUILD_TYPE=${build_type})"
