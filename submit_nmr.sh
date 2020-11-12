#!/bin/bash -l
#SBATCH -p normal
#SBATCH -J test
#SBATCH -n 16
#SBATCH --time=13-0:00:00
#SBATCH --mem-per-cpu=3300

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"


WORKDIR=/scratch/yanfeig/$SLURM_JOB_ID
export WORKDIR
mkdir -p $WORKDIR

export PATH=/opt/g16/:/opt/gv:$PATH
g16root=/opt export g16root
. /home/yanfeig/g16.profile

source /home/yanfeig/.bashrc
eval "$(conda shell.bash hook)"
conda activate /home/yanfeig/miniconda3/envs/nmr

python main.py
