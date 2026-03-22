module load nvhpc/25.3
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi/bin:$PATH
: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --cudamath-lib=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/lib64/ --use-mpi=yes --mpi-inc=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi/include/ --mpi-lib=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi/lib/ --use-nccl=yes --nccl-inc=$NVHPC_ROOT/comm_libs/nccl/include --nccl-lib=$NVHPC_ROOT/comm_libs/nccl/lib --use-cudss=yes --use-lua=yes;

