#!/bin/bash
##Resource Request
#SBATCH --job-name USER_JOB_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:L40S:1
#SBATCH --mem=USER_MEM
#SBATCH --time=USER_TIME
#SBATCH --partition=gpu

##Load the CUDA module
module load cuda
module load conda

conda init
source activate base
conda activate USER_CONDA_ENV

cellbender remove-background \
        --cuda \
        --input USER_RAW_MATRIX_FILE \
        # --expected-cells USER_EXPECTED_CELL_COUNT \
        --output USER_OUTPUT_MATRIX_FILE
