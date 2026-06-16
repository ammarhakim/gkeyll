module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.30
module load cudatoolkit/12.9
module load nccl/2.18.3-cu12
module load cray-libsci/25.09.0

: "${PREFIX:=/pscratch/sd/m/mana/gkeyll/code/gkeyll_gpu1/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-lib-name=sci_gnu --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib --cudamath-lib=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/math_libs/12.9/lib64 --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-nccl=yes --nccl-inc=$NCCL_DIR/include --nccl-lib=$NCCL_DIR/lib --use-lua=yes --use-cudss=yes;
