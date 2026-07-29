#!/bin/bash
#SBATCH --job-name=shearGPU
#SBATCH --output=/scratch/%u/matlab_logs/job_%A_%a.out
#SBATCH --error=/scratch/%u/matlab_logs/job_%A_%a.out
#SBATCH --array=1-1000%100
#SBATCH --partition=capacity
#SBATCH --gres=gpu:l4:1
#SBATCH --cpus-per-task=16 # This is cluster policy 16 CPU per GPU
#SBATCH --mem=3G
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=c.kawamura@uva.nl
#SBATCH --exclude=hipster-cn011

module load matlab/r2024a
# Force matlab to only use 1 thread per library
# Usually matlab will try to parallellize CPU math across all 16 CPU's 
# but since GPU is2026-06-27oing all the heavy lifting, this just wastes resource
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /home/ckawamu/repos/GranMA

CMDFILE="${CMDFILE:-hpc_scripts/commandsShearGPU.txt}"
IDX=$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))
CMD=$(sed -n "${IDX}p" "$CMDFILE")

echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID idx=$IDX] CMD: $CMD"

if [ -z "$CMD" ]; then
  echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID] ERROR: Empty command, exiting."
  exit 1
fi

eval "$CMD"


