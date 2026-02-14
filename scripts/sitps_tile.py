#!/usr/bin/env python3

"""Thin wrapper for the sitps_tile C++ binary."""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys
from typing import List


def parse_positive_int(raw: str) -> int:
    try:
        value = int(raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid integer value: {raw}") from exc
    if value <= 0:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {raw}")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tile/replicate SplitIndexTPS data.")
    parser.add_argument("--input-dir", required=True, help="Input SplitIndexTPS directory")
    parser.add_argument("--output-dir", required=True, help="Output SplitIndexTPS directory")
    parser.add_argument("--target-ly", required=True, type=parse_positive_int, help="Target lattice Ly")
    parser.add_argument("--target-lx", required=True, type=parse_positive_int, help="Target lattice Lx")
    parser.add_argument("--unit-ly", type=parse_positive_int, help="Optional PBC unit cell Ly")
    parser.add_argument("--unit-lx", type=parse_positive_int, help="Optional PBC unit cell Lx")
    parser.add_argument("--physics-json", help="Optional physics json for cross-check")
    return parser.parse_args()


def resolve_binary() -> str:
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    candidates = [
        repo_root / "build" / "sitps_tile",
        repo_root / "build_llvm_sdk" / "sitps_tile",
        repo_root / "sitps_tile",
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_mode & 0o111:
            return str(path)

    in_path = shutil.which("sitps_tile")
    if in_path:
        return in_path

    raise FileNotFoundError(
        "sitps_tile binary not found. Build target first, e.g. `cmake --build build --target sitps_tile`."
    )


def build_command(binary: str, args: argparse.Namespace) -> List[str]:
    cmd = [
        binary,
        "--input-dir",
        args.input_dir,
        "--output-dir",
        args.output_dir,
        "--target-ly",
        str(args.target_ly),
        "--target-lx",
        str(args.target_lx),
    ]
    if args.unit_ly is not None:
        cmd.extend(["--unit-ly", str(args.unit_ly)])
    if args.unit_lx is not None:
        cmd.extend(["--unit-lx", str(args.unit_lx)])
    if args.physics_json:
        cmd.extend(["--physics-json", args.physics_json])
    return cmd


def main() -> int:
    args = parse_args()
    try:
        binary = resolve_binary()
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 127

    cmd = build_command(binary, args)
    completed = subprocess.run(cmd, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
