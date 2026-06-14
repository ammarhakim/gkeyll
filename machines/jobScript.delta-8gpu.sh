#!/bin/bash
#SBATCH --job-name=fh-8gpu-test
#SBATCH --output=fh-8gpu-test.o%j
#SBATCH --error=fh-8gpu-test.e%j
#SBATCH --partition=gpuA100x4
#SBATCH --exclusive
#SBATCH --mem=0		
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --account=
#SBATCH --no-requeue
#SBATCH -t 00:05:00
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

module load nvhpc/25.3

srun -n 8 /projects/betv/gkylsoft/gkeyll/bin/gkeyll vlasov_par_firehose_beta3_1x3v_p2_L30_mu100_NX256_8GPUtest.lua



