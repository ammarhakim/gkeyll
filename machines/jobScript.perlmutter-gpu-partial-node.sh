#!/bin/bash -l

#.This jobscript is for partially full nodes (not using all the GPUs/node).

#.Declare a name for this job, preferably 16 or fewer characters.
#SBATCH -J <Job Name>

#.Enter the account to charge.
#SBATCH -A <Account Number>

#.Specify a queue.
#SBATCH -q regular

#.Number of nodes to request (Perlmutter 4 GPUs per node).
#SBATCH -N 2

#.Number of MPI processes per node (and we match 1 process to 1 GPU).
#SBATCH --ntasks-per-node 3

#.CPUs per MPI process, see Perlmutter docs.
#SBATCH --cpus-per-task 42

#.Specify GPU needs:
#SBATCH --constraint gpu
#SBATCH --gpus-per-node 4

#.Request wall time
#SBATCH -t 00:30:00

#.Mail is sent to you when the job starts and when it terminates or aborts.
#SBATCH --mail-user=<your email>
#SBATCH --mail-type=END,FAIL,REQUEUE

#.Load modules (this must match those in the machines/configure script).
module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/9.0.1
module load cudatoolkit/13.0
module load nccl/2.29.2-cu13
module load cray-libsci/25.09.0

#.Disable CUDA-ware MPI, since it causes problems on Perlmutter and we use NCCL alone.
export MPICH_GPU_SUPPORT_ENABLED=0

#.On Perlmutter some jobs get warnings about DVS_MAXNODES (used in file stripping).
#.We set it to 24 for now, but really this depends on the amount/size of I/O being performed.
#.See online NERSC docs and the intro_mpi man page.
export DVS_MAXNODES=24_
export MPICH_MPIIO_DVS_MAXNODES=24

# Safely route GPUDirect RDMA over the Host Bridge
export NCCL_NET_GDR_LEVEL=PHB

# Tell NCCL to use the Libfabric plugin
export NCCL_NET="AWS Libfabric"
export NCCL_CROSS_NIC=1

# Disable host registration and eager messages to prevent Slingshot 11 hangs
export FI_CXI_DISABLE_HOST_REGISTER=1
export FI_CXI_RDZV_GET_MIN=0
export FI_CXI_RDZV_THRESHOLD=0
export FI_CXI_RDZV_EAGER_SIZE=0

#.Run the rt_gk_sheath_2x2v_p1 executable using 1 GPU along x (-c 1) and 6
#.GPUs along the field line (-d 6). See './rt_gk_sheath_2x2v_p1 -h' for
#.more details/options on decomposition. It also assumes the executable is
#.in the present directory. If it isn't, change `./` to point to the
#.directory containing the executable.

echo "srun -u -N 2 -n 6 --ntasks-per-node 3 --gpus-per-node 4 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 6"
srun -u -N 2 -n 6 --ntasks-per-node 3 --gpus-per-node 4 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 6




