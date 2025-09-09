# Loading modules, pitagora currently only supports cuda 12.6

module load cuda/12.6
module load openmpi/4.1.6--gcc--12.3.0
module load nccl/2.22.3-1--gcc--12.3.0
module load cmake/3.27.9

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=90 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-mpi=yes --mpi-inc=$OPENMPI_HOME/include --mpi-lib=$OPENMPI_HOME/lib --use-nccl=yes --nccl-inc=$NCCL_HOME/include --nccl-lib=$NCCL_HOME/lib --use-lua=yes --use-cudss=yes;
