#. This is the loaded modules when the compilation and test succeded on Perlmutter GPU (ACDH 02/2025)
#.  1) craype-x86-milan                                5) gcc-native/13.2          9) sqs/2.0               13) cray-libsci/24.07.0             17) nccl/2.18.3-cu12 (g)
#.  2) libfabric/1.20.1                                6) perftools-base/24.07.0  10) darshan/default       14) PrgEnv-gnu/8.5.0      (prgenv)  18) cudatoolkit/12.4 (g)
#.  3) craype-network-ofi                              7) cpe/24.07               11) craype/2.7.32    (c)  15) craype-accel-nvidia80 (cpe)
#.  4) xpmem/2.9.6-1.1_20240510205610__g087dc11fc19d   8) gpu/1.0                 12) cray-dsmml/0.3.0      16) cray-mpich/8.1.28     (mpi)
module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.4
module load nccl/2.18.3-cu12

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-nccl=yes --nccl-inc=$NCCL_DIR/include --nccl-lib=$NCCL_DIR/lib --use-lua=yes --use-cudss=yes;
