module load nvhpc/25.3
: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --cudamath-lib=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/lib64/ --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include/ --mpi-lib=$CRAY_MPICH_DIR/lib/ --use-nccl=yes --nccl-inc=$NVHPC_ROOT/comm_libs/nccl/include --nccl-lib=$NVHPC_ROOT/comm_libs/nccl/lib --use-cudss=yes --use-lua=yes;

