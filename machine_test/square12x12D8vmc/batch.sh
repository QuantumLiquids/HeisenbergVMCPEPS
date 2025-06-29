#!/bin/bash

#SBATCH --job-name square12x12D8
#SBATCH --output=slurm-%j-%x.log
#SBATCH --get-user-env
#SBATCH --partition=normal
##SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
##SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --ntasks 100
#SBATCH --mem=800G  
##SBATCH --mem-per-cpu=8G
##SBATCH -C Intel
##SBATCH --nodes=4
##SBATCH --exclusive
##SBATCH --nodelist=c3
#SBATCH --dependency=46782

host=$(hostname)
echo "Hostname: $host"


source /opt/intel/oneapi/mkl/latest/env/vars.sh
source /opt/intel/oneapi/compiler/latest/env/vars.sh
export LIBRARY_PATH=$LD_LIBRARY_PATH

currTime=$(date +"%Y-%m-%d %T")
echo $currTime


rundir=$(pwd)
rootdir="../../"
executabledir="${rootdir}/build/"

#unset I_MPI_PMI_LIBRARY 
#ulimit -c unlimited
#ulimit -a
#export I_MPI_DEBUG=5
#${executabledir}/simple_update simple_update_params.json
mpirun -np ${SLURM_NTASKS} ${executabledir}/square_vmc_update vmc_params.json >> ${rundir}/${SLURM_JOB_NAME}.log

