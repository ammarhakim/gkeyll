#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks 4
#SBATCH --gpus 4 
#SBATCH --cpus-per-task=9
#SBATCH --partition=gpu_a100
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[mail address]

## Make sure the following lines (with # removed) are in ~/.bashrc:
# module load 2023													
# module load CUDA/12.1.1
# module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-12.1.1
# module load NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1
##

cd $HOME/gkylzero

srun -u -n 4 --gpus 4 ./cuda-build/rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 4