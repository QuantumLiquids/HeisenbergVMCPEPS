---
name: cluster-info
description: "Reference for susphy HPC cluster connection, hardware, and environment.
  Use when SSH'ing, submitting jobs, or diagnosing cluster issues."
---

# Susphy HPC Cluster Reference

## SSH Connection

Always use the local SSH alias:

```bash
ssh susphy
```

### Local SSH config (~/.ssh/config)

```
Host susphy susphy-login3
  HostName 10.20.26.130
  Port 22
  User wanghx
  ProxyCommand nc -X 5 -x 127.0.0.1:1080 %h %p
```

Requires a SOCKS5 proxy on localhost:1080 (ClashX or similar).

### Login node

Login lands on `login3.sharehpc.cn`. The alternate login4 is at 10.20.8.67.

## Home Directory

```
/share/home/wanghx/
```

Project repo: `~/HeisenbergVMCPEPS/`

## Hardware (256G56c partition)

- **56 cores per node** (2x Intel Xeon), 256 GB RAM
- All production runs use this partition
- The 24c partitions have ABI incompatibility with login-node builds

## Multi-Node Runs

- 2 nodes = 112 ranks, 512 GB total
- `MC_total_samples` must be divisible by `ntasks`
- 16x16 D=8 requires 2 nodes (single-node OOM at ~220 GB)

## Environment

Slurm scripts must source the environment with nounset temporarily disabled:

```bash
set +u; set +x
source /share/home/wanghx/myenv.sh
set -x; set -u
```

### Compilers

- **GCC 12.1** (required for C++20): `/share/software/gcc/12.1.0/bin/g++`
- GCC 8.5 (system default): lacks `[[unlikely]]` support, do not use
- hptt must be GCC-compiled (Intel-compiled hptt causes `_intel_fast_memcpy` link errors)

### Build on cluster

Install TensorToolkit, UltraDMRG, and PEPS into a common prefix in dependency
order. Build PEPS with `-DQLPEPS_USE_SCALAPACK=ON` when distributed MinSR is
required. The example below assumes that prefix is `$HOME/.local`.

```bash
cd ~/HeisenbergVMCPEPS/build
cmake .. \
  -DCMAKE_CXX_COMPILER=/share/software/gcc/12.1.0/bin/g++ \
  -DCMAKE_C_COMPILER=/share/software/gcc/12.1.0/bin/gcc \
  -DCMAKE_PREFIX_PATH=$HOME/.local
make -j4
```

## SCP File Transfer

```bash
# Upload to cluster
scp local_file susphy:~/HeisenbergVMCPEPS/path/

# Download from cluster
scp susphy:~/HeisenbergVMCPEPS/path/file .

# For directories: create remote dir first, then copy files
ssh susphy "mkdir -p ~/HeisenbergVMCPEPS/target/"
scp local_dir/* susphy:~/HeisenbergVMCPEPS/target/
```

Note: `scp -r` can fail if the remote target doesn't exist. Always `mkdir -p` first.

## Slurm Quick Reference

```bash
# Check running jobs
ssh susphy "squeue -u wanghx"

# Job details (memory, runtime)
ssh susphy "sacct -j JOBID --format=JobID,MaxRSS,Elapsed,State"

# Cancel a job
ssh susphy "scancel JOBID"

# Submit with dependency
ssh susphy "cd ~/HeisenbergVMCPEPS/run/... && sbatch --dependency=afterok:PREV_ID script.slurm"
```

## Known Issues

- **Shebang in heredoc**: Creating slurm scripts via `cat << 'EOF'` over SSH
  can produce `#\!/bin/bash` instead of `#!/bin/bash`. Always verify or fix with
  `sed -i '1s|.*|#!/bin/bash|'`.
- **Git pull conflicts**: SCP'd files may conflict with incoming commits.
  Remove untracked copies before `git pull`.
- **Proxy drops**: If `ssh susphy` fails with "Connection closed by 127.0.0.1
  port 1080/7897", the SOCKS proxy tunnel is down. User must restart it.
