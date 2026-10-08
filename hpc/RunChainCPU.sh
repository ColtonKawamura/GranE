#!/bin/bash
# CPU-based submission script for packRepeatTile jobs (repack).
# Modeled on RunChainGPU.sh but without GPU requirements so it can run
# concurrently with GPU sim jobs on the capacity partition.
#SBATCH --job-name=repack
#SBATCH --output=/scratch/%u/matlab_logs/job_%A_%a.out
#SBATCH --error=/scratch/%u/matlab_logs/job_%A_%a.out
#SBATCH --array=1-3%3
#SBATCH --partition=capacity
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=c.kawamura@uva.nl
#SBATCH --exclude=hipster-cn011

module load matlab/r2024a
# Force matlab to only use 1 thread per library
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /home/ckawamu/repos/GranE

CMDFILE="${CMDFILE:-hpc/commandsRepPack.txt}"
IDX=$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))
echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID] PWD at job start: $PWD"
echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID] SLURM_SUBMIT_DIR: $SLURM_SUBMIT_DIR"
CMD=$(sed -n "${IDX}p" "$CMDFILE")

echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID idx=$IDX] CMD: $CMD"

if [ -z "$CMD" ]; then
  echo "[job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID] ERROR: Empty command, exiting."
  exit 1
fi

eval "$CMD"
